#!/usr/bin/env python

# This file is part of Cantera. See License.txt in the top-level directory or
# at https://www.cantera.org/license.txt for license and copyright information.

"""
cti2yaml.py: Convert legacy CTI input files to YAML

Usage:
    python ck2yaml.py mech.cti [out.yaml]

This will produce the output file 'mech.yaml' if an output file name is not
specified.
"""

from __future__ import print_function

import sys
import os
import re
from collections import OrderedDict
import numpy as np

try:
    import ruamel_yaml as yaml
except ImportError:
    from ruamel import yaml

# Python 2/3 compatibility
try:
  basestring
except NameError:
  basestring = str


def _printerr(*args):
    # All debug and error output should go to stderr
    print(*args, file=sys.stderr)


class InputError(Exception):
    """
    Exception raised if an error is encountered while parsing the input file.
    """
    def __init__(self, msg, *args):
        if args:
            msg = msg.format(*args)
        super(InputError, self).__init__(msg)

BlockMap = yaml.comments.CommentedMap

def FlowMap(*args, **kwargs):
    m = yaml.comments.CommentedMap(*args, **kwargs)
    m.fa.set_flow_style()
    return m

def FlowList(*args, **kwargs):
    lst = yaml.comments.CommentedSeq(*args, **kwargs)
    lst.fa.set_flow_style()
    return lst

def represent_float(self, data):
    # type: (Any) -> Any
    if data != data:
        value = u'.nan'
    elif data == self.inf_value:
        value = u'.inf'
    elif data == -self.inf_value:
        value = u'-.inf'
    else:
        if data == 0:
            value = u'0.0'
        elif 0.01 <= abs(data) < 10000:
            value = np.format_float_positional(data, trim='0')
        else:
            value = np.format_float_scientific(data, trim='0')

    return self.represent_scalar(u'tag:yaml.org,2002:float', value)

yaml.RoundTripRepresenter.add_representer(float, represent_float)


def applyUnits(value):
    if isinstance(value, (float, int)):
        return value
    else:
        units = value[1]
        units = re.sub(r'([A-Za-z])-([A-Za-z])', r'\1*\2', units)
        units = re.sub(r'([A-Za-z])([-\d])', r'\1^\2', units)
        if value[0] == 0:
            return '0.0 {}'.format(units)
        elif 0.01 <= abs(value[0]) < 10000:
            return '{} {}'.format(
                np.format_float_positional(value[0], trim='0'), units)
        else:
            return '{} {}'.format(
                np.format_float_scientific(value[0], trim='0'), units)


# map of old CTI/XML names to the new YAML names
_newNames = {
    'GasKinetics': 'gas',
    'Interface': 'surface',
    'Edge': 'edge',
    'Mix': 'mixture-averaged',
    'Multi': 'multicomponent'
}

# constants that can be used in .cti files
OneAtm = 1.01325e5
OneBar = 1.0e5
# Conversion from eV to J/kmol (electronCharge * Navrog)
eV = 9.64853364595687e7
# Electron Mass in kg
ElectronMass = 9.10938291e-31

# default units
_ulen = 'm'
_umol = 'kmol'
_umass = 'kg'
_utime = 's'
_ue = 'J/kmol'
_uenergy = 'J'
_upres = 'Pa'

# default std state pressure
_pref = OneAtm

_name = 'noname'

# these lists store top-level entries
_elements = []
_species = []
_speciesnames = []
_phases = []
_reactions = {'reactions': []}

# default for Motz & Wise correction
_motz_wise = None

def enable_motz_wise():
    global _motz_wise
    _motz_wise = True

def disable_motz_wise():
    global _motz_wise
    _motz_wise = False

def validate(species = 'yes', reactions = 'yes'):
    pass

def dataset(nm):
    "Set the dataset name. Invoke this to change the name of the XML file."
    global _name
    _name = nm

def standard_pressure(p0):
    """Set the default standard-state pressure."""
    global _pref
    _pref = p0

def units(length = '', quantity = '', mass = '', time = '',
          act_energy = '', energy = '', pressure = ''):
    """
    Set the default units.

    :param length:
        The default units for length. Default: ``'m'``
    :param mass:
        The default units for mass. Default: ``'kg'``
    :param quantity:
        The default units to specify number of molecules. Default: ``'kmol'``
    :param time:
        The default units for time. Default: ``'s'``
    :param energy:
        The default units for energies. Default: ``'J'``
    :param act_energy:
        The default units for activation energies. Default: ``'K'``
    :param pressure:
        The default units for pressure. Default: ``'Pa'``
    """
    global _ulen, _umol, _ue, _utime, _umass, _uenergy, _upres
    if length: _ulen = length
    if quantity: _umol = quantity
    if act_energy: _ue = act_energy
    if time: _utime = time
    if mass: _umass = mass
    if energy: _uenergy = energy
    if pressure: _upres = pressure


def getComposition(atoms):
    if isinstance(atoms, dict): return atoms
    a = atoms.replace(',',' ')
    toks = a.split()
    d = OrderedDict()
    for t in toks:
        b = t.split(':')
        try:
            d[b[0]] = int(b[1])
        except ValueError:
            d[b[0]] = float(b[1])
    return d


class element(object):
    """ An atomic element or isotope. """
    def __init__(self, symbol='', atomic_mass=0.01, atomic_number=None):
        """
        :param symbol:
            The symbol for the element or isotope.
        :param atomic_mass:
            The atomic mass in amu.
        """
        self._sym = symbol
        self._atw = atomic_mass
        self._num = atomic_number
        _elements.append(self)

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('symbol', node._sym),
                        ('atomic-weight', node._atw)])
        if node._num is not None:
            out['atomic-number'] = node._num
        return representer.represent_dict(out)


class species(object):
    """A constituent of a phase or interface."""

    def __init__(self, name, atoms='', note = '', thermo=None, transport=None,
                 charge=None, size=1.0, standardState=None):
        """
        :param name:
            The species name (or formula). The name may be arbitrarily long,
            although usually a relatively short, abbreviated name is most
            convenient. Required parameter.
        :param atoms:
            The atomic composition, specified by a string containing
            space-delimited ``<element>:<atoms>`` pairs. The number of atoms may be
            either an integer or a floating-point number.
        :param note:
            A user-defined comment. Not evaluated by Cantera itself.
        :param thermo:
            The parameterization to use to compute the reference-state
            thermodynamic properties. This must be one of the entry types
            described in `Thermodynamic Property Models
            <https://cantera.org/science/science-species.html#sec-thermo-models>`_.
            To specify multiple parameterizations, each for a different temperature range,
            group them in parentheses.
        :param transport:
            An entry specifying parameters to compute this species'
            contribution to the transport properties. This must be one of the
            entry types described in `Species Transport Coefficients
            <https://cantera.org/science/science-species.html#species-transport-coefficients>`_,
            and must be consistent with the transport model of the phase into which
            the species is imported. To specify parameters for multiple
            transport models, group the entries in parentheses.
        :param size:
            The species "size". Currently used only for surface species,
            where it represents the number of sites occupied.
        :param charge:
            The charge, in multiples of :math:`|e|`. If not specified, the
            charge will be calculated from the number of "atoms" of element
            ``E``, which represents an electron.
        """
        self._name = name
        self._atoms = getComposition(atoms)
        if charge is not None and 'E' not in self._atoms:
            self._atoms['E'] = -charge
        self._size = size
        self._comment = note

        if isinstance(thermo, (list, tuple)):
            if isinstance(thermo[0], (NASA, NASA9, Shomate)):
                self._thermo = MultiPolyThermo(thermo)
        elif isinstance(thermo, (NASA, NASA9, Shomate)):
            self._thermo = MultiPolyThermo([thermo])
        elif thermo is not None:
            self._thermo = thermo
        else:
            self._thermo = const_cp()

        self._transport = transport
        self._standardState = standardState

        self.rk_pure = {}
        self.rk_binary = {}
        self.density = None

        _species.append(self)
        _speciesnames.append(name)

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('name', node._name),
                        ('composition', FlowMap(node._atoms.items()))])
        if node._size != 1:
            out['sites'] = node._size
        out['thermo'] = node._thermo
        if node.density:
            out['equation-of-state'] = {
                'model': 'constant-volume',
                'density': applyUnits(node.density)
            }

        if node.rk_pure:
            a = node.rk_pure['a']
            if isinstance(a, (tuple, list)):
                a = FlowList([applyUnits(ai) for ai in a])
            else:
                a = applyUnits(a)
            out['equation-of-state'] = {
                'model': 'Redlich-Kwong',
                'a': a,
                'b': applyUnits(node.rk_pure['b'])
            }

        if node.rk_binary:
            rkbin = BlockMap()
            for species, a in node.rk_binary.items():
                if isinstance(a, (tuple, list)):
                    rkbin[species] = FlowList([applyUnits(ai) for ai in a])
                else:
                    rkbin[species] = applyUnits(a)
            out['equation-of-state']['binary-a'] = rkbin

        if node._standardState:
            out['equation-of-state'] = {
                'model': 'constant-volume',
                'molar-volume': applyUnits(node._standardState._mv)
            }

        if node._transport:
            out['transport'] = node._transport
        if node._comment:
            out['note'] = node._comment
        return representer.represent_dict(out)


class thermo(object):
    """Base class for species thermodynamic properties."""
    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap()
        node.get_yaml(out)
        return representer.represent_dict(out)

    def get_yaml(self, out):
        out['model'] = self.model
        pref = self.pref or _pref
        if pref != OneAtm:
            out['reference-pressure'] = pref


class NASA(thermo):
    """The 7-coefficient NASA polynomial parameterization."""
    def __init__(self, Trange=(0.0, 0.0), coeffs=[], p0=None):
        r"""
        :param Trange:
            The temperature range over which the parameterization is valid.
            This must be entered as a sequence of two temperature values.
            Required.
        :param coeffs:
            List of seven coefficients :math:`(a_0, \ldots , a_6)`
        :param p0:
            The reference-state pressure, usually 1 atm or 1 bar. If omitted,
            the default value is used, which is set by the ``standard_pressure``
            directive.
        """
        self.model = 'NASA7'
        self._t = Trange
        self.pref = p0
        if len(coeffs) != 7:
            raise InputError('NASA coefficient list must have length = 7')
        self._coeffs = coeffs


class NASA9(thermo):
    """NASA9 polynomial parameterization for a single temperature region."""

    def __init__(self, Trange=(0.0, 0.0), coeffs=[], p0=None):
        r"""
        :param Trange:
            The temperature range over which the parameterization is valid.
            This must be entered as a sequence of two temperature values.
            Required.
        :param coeffs:
            List of nine coefficients :math:`(a_0, \ldots , a_8)`
        :param p0:
            The reference-state pressure, usually 1 atm or 1 bar. If omitted,
            the default value is used, which is set by the ``standard_pressure``
            directive.
        """
        self.model = 'NASA9'
        self._t = Trange         # Range of the polynomial representation
        self.pref = p0          # Reference pressure
        if len(coeffs) != 9:
            raise InputError('NASA9 coefficient list must have length = 9')
        self._coeffs = coeffs


class MultiPolyThermo(thermo):
    def __init__(self, regions):
        regions = sorted(regions, key=lambda r: r._t[0])
        self.pref = regions[0].pref
        self.Tranges = [regions[0]._t[0]]
        self.model = regions[0].model
        self.data = []
        for r in regions:
            self.Tranges.append(r._t[1])
            self.data.append(r._coeffs)

    def get_yaml(self, out):
        super(MultiPolyThermo, self).get_yaml(out)
        out['temperature-ranges'] = FlowList(self.Tranges)
        out['data'] = [FlowList(coeffs) for coeffs in self.data]


class Shomate(thermo):
    """Shomate polynomial parameterization."""
    def __init__(self, Trange=(0.0, 0.0), coeffs=[], p0=None):
        r"""
        :param Trange:
            The temperature range over which the parameterization is valid.
            This must be entered as a sequence of two temperature values.
            Required input.
        :param coeffs:
            Sequence of seven coefficients :math:`(A, \ldots ,G)`
        :param p0:
            The reference-state pressure, usually 1 atm or 1 bar. If omitted,
            the default value set by the ``standard_pressure`` directive is used.
        """
        self.model = 'Shomate'
        self._t = Trange
        self.pref = p0
        if len(coeffs) != 7:
            raise InputError('Shomate coefficient list must have length = 7')
        self._coeffs = coeffs


class const_cp(thermo):
    """Constant specific heat."""

    def __init__(self, t0=None, cp0=None, h0=None, s0=None, tmax=None,
                 tmin=None):
        """
        :param t0:
            Temperature parameter T0. Default: 298.15 K.
        :param cp0:
            Reference-state molar heat capacity (constant). Default: 0.0.
        :param h0:
            Reference-state molar enthalpy at temperature T0. Default: 0.0.
        :param s0:
            Reference-state molar entropy at temperature T0. Default: 0.0.
        """
        self.model = 'constant-cp'
        self.pref = None
        self.t0 = t0
        self.h0 = h0
        self.s0 = s0
        self.cp0 = cp0

    def get_yaml(self, out):
        super(const_cp, self).get_yaml(out)
        if self.t0 is not None:
            out['T0'] = applyUnits(self.t0)
        if self.h0 is not None:
            out['h0'] = applyUnits(self.h0)
        if self.s0 is not None:
            out['s0'] = applyUnits(self.s0)
        if self.cp0 is not None:
            out['cp0'] = applyUnits(self.cp0)


class gas_transport(object):
    """
    Species-specific Transport coefficients for gas-phase transport models.
    """
    def __init__(self, geom, diam, well_depth, dipole=0.0, polar=0.0,
                 rot_relax=0.0, acentric_factor=None, disp_coeff=0.0,
                 quad_polar=0.0):
        """
        :param geom:
            A string specifying the molecular geometry. One of ``atom``,
            ``linear``, or ``nonlinear``. Required.
        :param diam:
            The Lennard-Jones collision diameter in Angstroms. Required.
        :param well_depth:
            The Lennard-Jones well depth in Kelvin. Required.
        :param dipole:
            The permanent dipole moment in Debye. Default: 0.0
        :param polar:
            The polarizability in A^3. Default: 0.0
        :param rot_relax:
            The rotational relaxation collision number at 298 K. Dimensionless.
            Default: 0.0
        :param w_ac:
            Pitzer's acentric factor.  Dimensionless.
            Default: 0.0
        :param disp_coeff:
            The dispersion coefficient in A^5
            Default: 0.0
        :param quad_polar:
            The quadrupole polarizability
            Default: 0.0
        """
        self._geom = geom
        self._diam = diam
        self._well_depth = well_depth
        self._dipole = dipole
        self._polar = polar
        self._rot_relax = rot_relax
        self._w_ac = acentric_factor
        self._disp_coeff = disp_coeff
        self._quad_polar = quad_polar

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap([('model', 'gas'),
                        ('geometry', node._geom),
                        ('diameter', node._diam),
                        ('well-depth', node._well_depth)])
        if node._dipole:
            out['dipole'] = node._dipole
        if node._polar:
            out['polarizability'] = node._polar
        if node._rot_relax:
            out['rotational-relaxation'] = node._rot_relax
        if node._w_ac:
            out['acentric-factor'] = node._w_ac
        if node._disp_coeff:
            out['dispersion-coefficient'] = node._disp_coeff
        if node._quad_polar:
            out['quadrupole-polarizability'] = node._quad_polar
        return representer.represent_dict(out)


class Arrhenius(object):
    def __init__(self, A=0.0, b=0.0, E=0.0, coverage=[]):
        """
        :param A:
            The pre-exponential coefficient. Required input. If entered without
            units, the units will be computed considering all factors that
            affect the units. The resulting units string is written to the CTML
            file individually for each reaction pre-exponential coefficient.
        :param b:
            The temperature exponent. Dimensionless. Default: 0.0.
        :param E:
            Activation energy. Default: 0.0.
        :param coverage: For a single coverage dependency, a list with four
            elements: the species name followed by the three coverage
            parameters. For multiple coverage dependencies, a list of lists
            containing the individual sets of coverage parameters. Only used for
            surface and edge reactions.
        """

        self._c = [A, b, E]

        if coverage:
            if isinstance(coverage[0], basestring):
                self._cov = [coverage]
            else:
                self._cov = coverage
            for cov in self._cov:
                if len(cov) != 4:
                    raise InputError("Incorrect number of coverage parameters")
        else:
            self._cov = None

    @classmethod
    def to_yaml(cls, representer, node):
        out = FlowMap([('A', applyUnits(node._c[0])),
                       ('b', applyUnits(node._c[1])),
                       ('Ea', applyUnits(node._c[2]))])
        return representer.represent_dict(out)


class stick(Arrhenius):
    def __init__(self, *args, **kwargs):
        """
        :param motz_wise: 'True' if the Motz & Wise correction should be used,
            'False' if not. If unspecified, use the mechanism default (set using
            the functions `enable_motz_wise` or `disable_motz_wise`).
        """
        self.motz_wise = kwargs.pop('motz_wise', None)
        Arrhenius.__init__(self, *args, **kwargs)


class reaction(object):
    """
    A homogeneous chemical reaction with pressure-independent rate coefficient
    and mass-action kinetics.
    """
    def __init__(self, equation, kf, id='', order='', options=[]):
        r"""
        :param equation:
            A string specifying the chemical equation.
        :param kf:
            The rate coefficient for the forward direction. If a sequence of
            three numbers is given, these will be interpreted as [A, b, E] in
            the modified Arrhenius function :math:`A T^b exp(-E/\hat{R}T)`.
        :param id:
            An optional identification string. If omitted, it defaults to a
            four-digit numeric string beginning with 0001 for the first
            reaction in the file.
        :param order:
            Override the default reaction orders implied by the reactant
            stoichiometric coefficients. Given as a string of key:value pairs,
            e.g., ``"CH4:0.25 O2:1.5"``.
        :param options: Processing options, as described in
            `Options <https://cantera.org/tutorials/cti/reactions.html#options>`_.
            May be one or more (as a list) of the
            following: ``'duplicate'``, ``'negative_A'``,`` 'negative_orders'``,
            ``'nonreactant_orders'``.
        """
        self._e = equation
        self._order = getComposition(order)
        self._num = len(_reactions['reactions']) + 1
        self._id = id
        self._options = [options] if isinstance(options, str) else options
        self._kf = Arrhenius(*kf) if isinstance(kf, (list, tuple)) else kf
        self._type = 'elementary'
        _reactions['reactions'].append(self)

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap()
        node.get_yaml(out)
        return representer.represent_dict(out)

    def get_yaml(self, out):
        out['equation'] = self._e
        out.yaml_add_eol_comment('Reaction {}'.format(self._num), 'equation')
        if self._type not in ('elementary', 'edge', 'surface'):
            out['type'] = self._type

        if self._type in ('elementary', 'three-body', 'edge', 'surface'):
            out['rate-constant'] = self._kf

        if 'duplicate' in self._options:
            out['duplicate'] = True
        if 'negative_A' in self._options:
            out['negative-A'] = True

        if self._order:
            out['orders'] = FlowMap(self._order.items())
        if 'negative_orders' in self._options:
            out['negative-orders'] = True
        if 'nonreactant_orders' in self._options:
            out['nonreactant-orders'] = True


class three_body_reaction(reaction):
    """
    A three-body reaction.
    """
    def __init__(self, equation, kf, efficiencies='', id='', options=[]):
        """
        :param equation:
            A string specifying the chemical equation. The reaction can be
            written in either the association or dissociation directions, and
            may be reversible or irreversible.
        :param kf:
            The rate coefficient for the forward direction. If a sequence of
            three numbers is given, these will be interpreted as [A, b, E] in
            the modified Arrhenius function.
        :param efficiencies:
            A string specifying the third-body collision efficiencies.
            The efficiencies for unspecified species are set to 1.0.
        :param id:
            An optional identification string. If omitted, it defaults to a
            four-digit numeric string beginning with 0001 for the first
            reaction in the file.
        :param options: Processing options, as described in
            `Options <https://cantera.org/tutorials/cti/reactions.html#options>`_.
        """
        reaction.__init__(self, equation, kf, id, '', options)
        self._type = 'three-body'
        self._eff = getComposition(efficiencies)

    def get_yaml(self, out):
        super(three_body_reaction, self).get_yaml(out)
        if self._eff:
            out['efficiencies'] = FlowMap(self._eff)


class pdep_reaction(reaction):
    """ Base class for falloff_reaction and chemically_activated_reaction """
    def __init__(self, equation, klow, khigh, efficiencies, falloff, id, options):
        super(pdep_reaction, self).__init__(equation, None, id, '', options)
        self._klow = Arrhenius(*klow) if isinstance(klow, (list, tuple)) else klow
        self._khigh = Arrhenius(*khigh) if isinstance(khigh, (list, tuple)) else khigh
        self._falloff = falloff
        self._eff = getComposition(efficiencies)

    def get_yaml(self, out):
        super(pdep_reaction, self).get_yaml(out)

        out['low-P-rate-constant'] = self._klow
        out['high-P-rate-constant'] = self._khigh

        if self._falloff:
            self._falloff.get_yaml(out)

        if self._eff:
            out['efficiencies'] = FlowMap(self._eff)


class falloff_reaction(pdep_reaction):
    """ A gas-phase falloff reaction. """
    def __init__(self, equation, kf0, kf, efficiencies='', falloff=None, id='',
                 options=[]):
        """
        :param equation:
            A string specifying the chemical equation.
        :param kf:
            The rate coefficient for the forward direction in the high-pressure
            limit. If a sequence of three numbers is given, these will be
            interpreted as [A, b, E] in the modified Arrhenius function.
        :param kf0:
            The rate coefficient for the forward direction in the low-pressure
            limit. If a sequence of three numbers is given, these will be
            interpreted as [A, b, E] in the modified Arrhenius function.
        :param efficiencies:
            A string specifying the third-body collision efficiencies. The
            efficiency for unspecified species is set to 1.0.
        :param falloff:
            An embedded entry specifying a falloff function. If omitted, a
            unity falloff function (Lindemann form) will be used.
        :param id:
            An optional identification string. If omitted, it defaults to a
            four-digit numeric string beginning with 0001 for the first
            reaction in the file.
        :param options:
            Processing options, as described in
            `Options <https://cantera.org/tutorials/cti/reactions.html#options>`_.
        """
        super(falloff_reaction, self).__init__(equation, kf0, kf, efficiencies,
            falloff, id, options)
        self._type = 'falloff'


class chemically_activated_reaction(pdep_reaction):
    """ A gas-phase, chemically activated reaction. """

    def __init__(self, equation, kLow, kHigh,
                 efficiencies='', falloff=None, id='', options=[]):
        """
        :param equation:
            A string specifying the chemical equation.
        :param kLow:
            The rate coefficient for the forward direction in the low-pressure
            limit. If a sequence of three numbers is given, these will be
            interpreted as [A, b, E] in the modified Arrhenius function.
        :param kHigh:
            The rate coefficient for the forward direction in the high-pressure
            limit. If a sequence of three numbers is given, these will be
            interpreted as [A, b, E] in the modified Arrhenius function.
        :param efficiencies:
            A string specifying the third-body collision efficiencies. The
            efficiency for unspecified species is set to 1.0.
        :param falloff:
            An embedded entry specifying a falloff function. If omitted, a
            unity falloff function (Lindemann form) will be used.
        :param id:
            An optional identification string. If omitted, it defaults to a
            four-digit numeric string beginning with 0001 for the first
            reaction in the file.
        :param options:
            Processing options, as described in
            `Options <https://cantera.org/tutorials/cti/reactions.html#options>`_.
        """
        super(chemically_activated_reaction, self).__init__(equation, kLow,
            kHigh, efficiencies, falloff, id, options)
        self._type = 'chemically-activated'


class pdep_arrhenius(reaction):
    """
    Pressure-dependent rate calculated by interpolating between Arrhenius
    expressions at different pressures.

    :param equation:
        A string specifying the chemical equation.
    :param args:
        Each additional argument is a sequence of four elements specifying the
        pressure and the Arrhenius parameters at that pressure.
    """
    def __init__(self, equation, *args, **kwargs):
        reaction.__init__(self, equation, None, **kwargs)
        self.arrhenius = args
        self._type = 'pressure-dependent-Arrhenius'

    def get_yaml(self, out):
        super(pdep_arrhenius, self).get_yaml(out)
        rates = []
        for p, A, b, Ea in self.arrhenius:
            rates.append(FlowMap([('P', applyUnits(p)),
                                 ('A', applyUnits(A)),
                                 ('b', applyUnits(b)),
                                 ('Ea', applyUnits(Ea))]))
        out['rate-constants'] = rates


class chebyshev_reaction(reaction):
    """
    Pressure-dependent rate calculated in terms of a bivariate Chebyshev
    polynomial.

    :param equation:
        A string specifying the chemical equation.
    :param Tmin:
        The minimum temperature at which the rate expression is defined
    :param Tmax:
        the maximum temperature at which the rate expression is defined
    :param Pmin:
        The minimum pressure at which the rate expression is defined
    :param Pmax:
        The maximum pressure at which the rate expression is defined
    :param coeffs:
        A 2D array of the coefficients defining the rate expression. For a
        polynomial with M points in temperature and N points in pressure, this
        should be a list of M lists each with N elements.
    """
    def __init__(self, equation, Tmin=300.0, Tmax=2500.0, Pmin=(0.001, 'atm'),
                 Pmax=(100.0, 'atm'), coeffs=[[]], **kwargs):
        reaction.__init__(self, equation, None, **kwargs)
        self._type = 'Chebyshev'
        self.Pmin = Pmin
        self.Pmax = Pmax
        self.Tmin = Tmin
        self.Tmax = Tmax
        self.coeffs = coeffs

    def get_yaml(self, out):
        super(chebyshev_reaction, self).get_yaml(out)
        out['temperature-range'] = FlowList([applyUnits(self.Tmin),
                                             applyUnits(self.Tmax)])
        out['pressure-range'] = FlowList([applyUnits(self.Pmin),
                                          applyUnits(self.Pmax)])
        out['data'] = [FlowList(line) for line in self.coeffs]


class surface_reaction(reaction):
    """
    A heterogeneous chemical reaction with pressure-independent rate
    coefficient and mass-action kinetics.
    """
    def __init__(self, equation, kf, id='', order='', beta=None, options=[],
                 rate_coeff_type=''):
        """
        :param equation:
            A string specifying the chemical equation.
        :param kf:
            The rate coefficient for the forward direction. If a sequence of
            three numbers is given, these will be interpreted as [A, b, E] in
            the modified Arrhenius function.
        :param id:
            An optional identification string. If omitted, it defaults to a
            four-digit numeric string beginning with 0001 for the first
            reaction in the file.
        :param beta:
            Charge transfer coefficient: A number between 0 and 1 which, for a
            charge transfer reaction, determines how much of the electric
            potential difference between two phases is applied to the
            activation energy of the fwd reaction.  The remainder is applied to
            the reverse reaction.
        :param options:
            Processing options, as described in
            `Options <https://cantera.org/tutorials/cti/reactions.html#options>`_.
        """
        reaction.__init__(self, equation, kf, id, order, options)
        self._type = 'surface'
        self.sticking = isinstance(kf, stick)
        self._beta = beta
        self._rate_coeff_type = rate_coeff_type

    def get_yaml(self, out):
        super(surface_reaction, self).get_yaml(out)
        if self.sticking:
            del out['rate-constant']
            out.insert(1, 'sticking-coefficient', self._kf)
            if self._kf.motz_wise is not None:
                out['Motz-Wise'] = self._kf.motz_wise
        if self._rate_coeff_type == 'exchangecurrentdensity':
            out['exchange-current-density-formulation'] = True

        if self._kf._cov is not None:
            cov = {c[0]: FlowMap([('a', c[1]), ('m', c[2]), ('E', c[3])])
                   for c in self._kf._cov}
            out['coverage-dependencies'] = cov
        if self._beta is not None:
            out['beta'] = self._beta

        # todo: exchangecurrentdensity


class edge_reaction(surface_reaction):
    def __init__(self, equation, kf, id='', order='', beta=None, options=[],
                 rate_coeff_type=''):
        super(edge_reaction, self).__init__(equation, kf, id, order, beta,
              options, rate_coeff_type)
        self._type = 'edge'


class state(object):
    """
    An embedded entry that specifies the thermodynamic state of a phase
    or interface.
    """
    def __init__(self, temperature=None, pressure=None, mole_fractions=None,
                 mass_fractions=None, density=None, coverages=None,
                 solute_molalities=None):
        """
        :param temperature:
            The temperature.
        :param pressure:
            The pressure.
        :param density:
            The density. Cannot be specified if the phase is incompressible.
        :param mole_fractions:
            A string specifying the species mole fractions. Unspecified species
            are set to zero.
        :param mass_fractions:
            A string specifying the species mass fractions. Unspecified species
            are set to zero.
        :param coverages:
            A string specifying the species coverages. Unspecified species are
            set to zero. Can only be specified for interfaces.
        """
        self._t = temperature
        self._rho = density
        self._p = pressure
        self._x = mole_fractions
        self._y = mass_fractions
        self._c = coverages
        self._m = solute_molalities

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap()
        if node._t is not None:
            out['T'] = applyUnits(node._t)
        if node._p is not None:
            out['P'] = applyUnits(node._p)
        if node._rho is not None:
            out['density'] = applyUnits(node._rho)
        if node._x is not None:
            out['X'] = FlowMap(getComposition(node._x).items())
        if node._y is not None:
            out['Y'] = FlowMap(getComposition(node._y).items())
        if node._c is not None:
            out['coverages'] = FlowMap(getComposition(node._c).items())
        if node._m is not None:
            out['molalities'] = FlowMap(getComposition(node._m).items())
        return representer.represent_dict(out)


class phase(object):
    """Base class for phases of matter."""

    def __init__(self, name='', dim=3, elements='', species='', note='',
                 reactions='none', initial_state=None, options=[]):
        """
        :param name:
            A string to identify the phase. Must be unique among the phase
            names within the file.
        :param elements:
            The elements. A string of element symbols.
        :param species:
            The species. A string or sequence of strings in the format
            described in `Defining the Species
            <https://cantera.org/tutorials/cti/phases.html#defining-the-species>`_.
        :param note:
            A user-defined comment. Not evaluated by Cantera itself.
        :param reactions:
            The homogeneous reactions. If omitted, no reactions will be
            included. A string or sequence of strings in the format described
            in `Declaring the Reactions
            <https://cantera.org/tutorials/cti/phases.html#declaring-the-reactions>`_.
            This field is not allowed for ``stoichiometric_solid`` and
            ``stoichiometric_liquid`` entries.
        :param kinetics:
            The kinetics model. Optional; if omitted, the default model for the
            phase type will be used.
        :param transport:
            The transport property model. Optional. If omitted, transport
            property calculation will be disabled.
        :param initial_state:
            Initial thermodynamic state, specified with an embedded state entry.
        :param options:
            Special processing options. Optional.
        """

        self._name = name
        self._el = elements
        self._sp = []
        self._rx = []
        self._kin = None
        self._tr = None
        self._comment = note
        self._options = [options] if isinstance(options, str) else options

        #--------------------------------
        #        process species
        #--------------------------------

        # if a single string is entered, make it a list
        if isinstance(species, str):
            species = [species]

        # for each species string, check whether or not the species
        # are imported or defined locally. If imported, the string
        # contains a colon (:)
        for sp in species:
            foundColon = False
            allLocal = True
            for token in sp.split():
                if ':' in sp:
                    foundColon = True
                if token not in _speciesnames:
                    allLocal = False

            if foundColon and not allLocal:
                icolon = sp.find(':')
                datasrc = sp[:icolon].strip()
                spnames = sp[icolon+1:].strip()
                if spnames != 'all':
                    spnames = FlowList(spnames.split())
                self._sp.append((datasrc + '.yaml/species', spnames))

            else:
                spnames = sp
                self._sp.append(('species', FlowList(spnames.split())))

        if isinstance(reactions, str):
            reactions = [reactions]

        # for each reaction string, check whether or not the reactions
        # are imported or defined locally. If imported, the string
        # contains a colon (:)
        for r in reactions:
            icolon = r.find(':')
            if icolon > 0:
                datasrc = r[:icolon].strip() + '.yaml/reactions'
                rnum = r[icolon+1:].strip()
            else:
                datasrc = 'reactions'
                rnum = r.strip()
            if rnum == 'all' and 'skip_undeclared_species' in self._options:
                rnum = 'declared-species'
            if rnum != 'none':
                self._rx.append([datasrc, rnum])
            if rnum in ('all', 'declared-species', 'none'):
                continue
            if '*' in rnum:
                if datasrc != 'reactions':
                    _printerr("WARNING: Reaction id-pattern matching from remote"
                        " files not supported ({}: {})".format(datasrc, rnum))
            else:
                _printerr("WARNING: Reaction specification"
                          " '{}' not supported".format(rnum))

        self._initial = initial_state

        # add this phase to the global phase list
        _phases.append(self)

    @classmethod
    def to_yaml(cls, representer, node):
        out = BlockMap()
        node.get_yaml(out)
        return representer.represent_dict(out)

    def get_yaml(self, out):
        out['name'] = self._name
        out['thermo'] = self.thermo_model
        out['elements'] = FlowList(self._el.split())

        if len(self._sp) == 1 and self._sp[0][0] == 'species':
            # all local species
            out['species'] = self._sp[0][1]
        else:
            out['species'] = [BlockMap([(sp[0], sp[1])]) for sp in self._sp]

        if 'skip_undeclared_elements' in self._options:
            out['skip-undeclared-elements'] = True

        # Convert reaction pattern matching to use of multiple reaction sections
        for i in range(len(self._rx)):
            spec = self._rx[i][1]
            name = self._name + '-reactions'
            if '*' in spec and name not in _reactions:
                pattern = re.compile(spec.replace('*', '.*'))
                misses = []
                hits = []
                for reaction in _reactions['reactions']:
                    if pattern.match(reaction._id):
                        hits.append(reaction)
                    else:
                        misses.append(reaction)
                _reactions[name] = hits
                _reactions['reactions'] = misses
                self._rx[i] = [name, 'all']

        if self._kin and self._rx:
            out['kinetics'] = _newNames[self._kin]
            if len(self._rx) == 1 and self._rx[0][0] == 'reactions':
                out['reactions'] = self._rx[0][1]
            elif all(r[1] == 'all' for r in self._rx):
                out['reactions'] = FlowList(r[0] for r in self._rx)
            else:
                out['reactions'] = [BlockMap([(r[0], r[1])]) for r in self._rx]


        if self._tr:
            out['transport'] = _newNames[self._tr]

        if self._comment:
            out['note'] = self._comment

        if self._initial:
            out['state'] = self._initial


class ideal_gas(phase):
    """An ideal gas mixture."""
    def __init__(self, name='', elements='', species='', note='',
                 reactions='none', kinetics='GasKinetics', transport=None,
                 initial_state=None, options=[]):
        """
        The parameters correspond to those of :class:`.phase`, with the
        following modifications:

        :param kinetics:
            The kinetics model. Usually this field is omitted, in which case
            kinetics model GasKinetics, appropriate for reactions in ideal gas
            mixtures, is used.
        :param transport:
            The transport property model. One of the strings ``'none'``,
            ``'multi'``, or ``'mix'``. Default: ``'none'``.
        """

        phase.__init__(self, name, 3, elements, species, note, reactions,
                       initial_state, options)
        self._kin = kinetics
        self._tr = transport
        self.thermo_model = 'ideal-gas'


class stoichiometric_solid(phase):
    """
    A solid compound or pure element. Stoichiometric solid phases contain
    exactly one species, which always has unit activity. The solid is assumed
    to have constant density. Therefore the rates of reactions involving these
    phases do not contain any concentration terms for the (one) species in the
    phase, since the concentration is always the same."""
    def __init__(self, name='', elements='', species='', note='', density=None,
                 transport='None', initial_state=None, options=[]):
        """
        See :class:`.phase` for descriptions of the parameters.
        """

        phase.__init__(self, name, 3, elements, species, note, 'none',
                       initial_state, options)
        self.thermo_model = 'fixed-stoichiometry'
        self._dens = density
        if self._dens is None:
            raise InputError('density must be specified.')
        self._tr = None if transport == 'None' else transport

    def get_yaml(self, out):
        super(stoichiometric_solid, self).get_yaml(out)
        for section, names in self._sp:
            if section != 'species':
                _printerr("WARNING: Converting stoichiometric_solid species from"
                    " different input files ({}) is not supported.".format(section))
            else:
                species = [S for S in _species if S._name == names[0]][0]
                species.density = self._dens


class stoichiometric_liquid(stoichiometric_solid):
    """
    An incompressible stoichiometric liquid. Currently, there is no
    distinction between stoichiometric liquids and solids.
    """


class metal(phase):
    """A metal."""
    def __init__(self, name='', elements='', species='', note='', density=-1.0,
                 transport='None', initial_state=None, options=[]):

        phase.__init__(self, name, 3, elements, species, note, 'none',
                       initial_state, options)
        self.thermo_model = 'metal-electron'
        self._dens = density

    def get_yaml(self, out):
        super(metal, self).get_yaml(out)
        out['density'] = applyUnits(self._dens)


class incompressible_solid(phase):
    """An incompressible solid."""
    def __init__(self, name='', elements='', species='', note='', density=None,
                 transport='None', initial_state=None, options=[]):

        phase.__init__(self, name, 3, elements, species, note, 'none',
                       initial_state, options)
        self.thermo_model = 'constant-density'
        self._dens = density
        if self._dens is None:
            raise InputError('density must be specified.')

    def get_yaml(self, out):
        super(incompressible_solid, self).get_yaml(out)
        out['density'] = applyUnits(self._dens)


class liquid_vapor(phase):
    """
    A fluid with a complete liquid/vapor equation of state. This entry type
    selects one of a set of predefined fluids with built-in liquid/vapor
    equations of state. The substance_flag parameter selects the fluid. See
    liquidvapor.cti and liquidvapor.py for the usage of this entry type.
    """
    pure_fluids = {
        0: 'water',
        1: 'nitrogen',
        2: 'methane',
        3: 'hydrogen',
        4: 'oxygen',
        5: 'HFC134a',
        7: 'carbondioxide',
        8: 'heptane'
    }

    def __init__(self, name='', elements='', species='', note='',
                 substance_flag=0, initial_state=None, options=[]):

        phase.__init__(self, name, 3, elements, species, note, 'none',
                       initial_state, options)
        self.thermo_model = 'pure-fluid'
        self._subflag = substance_flag

    def get_yaml(self, out):
        super(liquid_vapor, self).get_yaml(out)
        if self._subflag in self.pure_fluids:
            out['pure-fluid-name'] = self.pure_fluids[self._subflag]
        else:
            raise InputError('liquid_vapor: unrecognized value "{}" for '
                '"substance_flag"', self._subflag)


class pureFluidParameters(object):
    def __init__(self, species=None, a_coeff=[], b_coeff=0):
        self._species = species
        self._acoeff = a_coeff
        self._bcoeff = b_coeff


class crossFluidParameters(object):
    def __init__(self, species=None, a_coeff=[], b_coeff=[]):
        self._species1, self._species2 = species.split(' ')
        self._acoeff = a_coeff
        self._bcoeff = b_coeff


class RedlichKwongMFTP(phase):
    """A multi-component fluid model for non-ideal gas fluids.
        """

    def __init__(self, name='', elements='', species='', note='',
                 reactions='none', kinetics='GasKinetics', initial_state=None,
                 activity_coefficients=None, transport='None', options=[]):

        phase.__init__(self,name, 3, elements, species, note, reactions,
                       initial_state,options)
        self.thermo_model = 'Redlich-Kwong'
        self._kin = kinetics
        self._tr = None if transport == 'None' else transport
        self._activityCoefficients = activity_coefficients

    def get_yaml(self, out):
        super(RedlichKwongMFTP, self).get_yaml(out)
        for section, names in self._sp:
            if section != 'species':
                _printerr("WARNING: Converting Redlich-Kwong species from"
                    " different input files ({}) is not supported.".format(section))

        spdict = {sp._name: sp for sp in _species}
        for params in self._activityCoefficients:
            if isinstance(params, pureFluidParameters):
                sp = spdict[params._species]
                sp.rk_pure = {'a': params._acoeff, 'b': params._bcoeff}
            elif isinstance(params, crossFluidParameters):
                sp1 = spdict[params._species1]
                sp1.rk_binary[params._species2] = params._acoeff
                sp2 = spdict[params._species2]
                sp2.rk_binary[params._species1] = params._acoeff


class constantIncompressible(object):
    """Constant molar volume."""
    def __init__(self, molarVolume=0.0):
        """
        :param molarVolume:
            Reference-state molar volume. Default: 0.0.
        """
        self._mv = molarVolume


class IdealSolidSolution(phase):
    """An IdealSolidSolution phase."""
    def __init__(self, name='', elements='', species='', note='',
                 transport='None', initial_state=None,
                 standard_concentration=None, options=[]):
        phase.__init__(self, name, 3, elements, species, note, 'none',
                       initial_state, options)
        self.thermo_model = 'ideal-solid-solution'
        self._stdconc = standard_concentration
        if self._stdconc is None:
            raise InputError('In phase {}: standard_concentration must be specified.', name)
        self._tr = None if transport == 'None' else transport

    def get_yaml(self, out):
        super(IdealSolidSolution, self).get_yaml(out)
        out['standard-concentration'] = self._stdconc.replace('_', '-')


class table(thermo):
    """User provided thermo table for BinarySolutionTabulatedThermo"""
    def __init__(self, moleFraction=([],''), enthalpy=([],''), entropy=([],'')):
        """
        :param moleFraction:
            The mole fraction of the tabulated species. Required parameter.
        :param enthalpy:
            The enthalpy of the tabulated species. Required parameter.
        :param entropy:
            The entropy of the tabulated species. Required parameter.
        """
        self.x = moleFraction
        self.h = enthalpy
        self.s = entropy


class BinarySolutionTabulatedThermo(IdealSolidSolution):
    """A BinarySolutionTabulatedThermo phase."""
    def __init__(self, name='', elements='', species='', note='',
                 transport='None', initial_state=None,
                 standard_concentration=None, tabulated_species=None,
                 tabulated_thermo=None, options=[]):
        IdealSolidSolution.__init__(self, name, elements, species, note,
                                    transport, initial_state,
                                    standard_concentration, options)
        self.thermo_model = 'binary-solution-tabulated'
        self._tabSpecies = tabulated_species
        self._tabThermo = tabulated_thermo
        self._stdconc = standard_concentration
        self._tr = None if transport == 'None' else transport
        if self._stdconc is None:
            raise InputError('In phase {}: standard_concentration must be specified.', name)
        if tabulated_species is None:
            raise InputError('In phase {}: tabulated_species must be specified.', name)
        if tabulated_thermo is None:
            raise InputError('In phase {}: Thermo data must be provided for the tabulated_species.', name)

    def get_yaml(self, out):
        super(BinarySolutionTabulatedThermo, self).get_yaml(out)
        out['tabulated-species'] = self._tabSpecies
        energy_units, quantity_units = self._tabThermo.h[1].split('/')
        tabThermo = BlockMap()
        if energy_units != _uenergy or quantity_units != _umol:
            tabThermo['units'] = FlowMap([('energy', energy_units),
                                          ('quantity', quantity_units)])
        tabThermo['mole-fractions'] = FlowList(self._tabThermo.x[0])
        tabThermo['enthalpy'] = FlowList(self._tabThermo.h[0])
        tabThermo['entropy'] = FlowList(self._tabThermo.s[0])
        out['tabulated-thermo'] = tabThermo

class ideal_interface(phase):
    """A chemically-reacting ideal surface solution of multiple species."""
    def __init__(self, name='', elements='', species='', note='',
                 reactions='none', site_density=0.0, phases=[],
                 kinetics='Interface', transport='None', initial_state=None,
                 options=[]):
        """
        The parameters correspond to those of :class:`.phase`, with the
        following modifications:

        :param reactions:
            The heterogeneous reactions at this interface. If omitted, no
            reactions will be included. A string or sequence of strings in the
            format described in `Declaring the Reactions
            <https://cantera.org/tutorials/cti/phases.html#declaring-the-reactions>`_.
        :param site_density:
            The number of adsorption sites per unit area.
        :param phases:
            A string listing the bulk phases that participate in reactions
            at this interface.
        """
        phase.__init__(self, name, 2, elements, species, note, reactions,
                       initial_state, options)
        self.thermo_model = 'surface'
        self._kin = kinetics
        self._tr = None if transport == 'None' else transport
        self._sitedens = site_density

    def get_yaml(self, out):
        super(ideal_interface, self).get_yaml(out)
        out['site-density'] = applyUnits(self._sitedens)


class edge(ideal_interface):
    """A 1D boundary between two surface phases."""
    def __init__(self, name='', elements='', species='', note='',
                 reactions='none', site_density=0.0, phases=[], kinetics='Edge',
                 transport='None', initial_state=None, options=[]):

        ideal_interface.__init__(self, name, elements, species, note, reactions,
            site_density, phases, kinetics, transport, initial_state, options)
        self.thermo_model = 'edge'


# Falloff parameterizations

class Troe(object):
    """The Troe falloff function."""
    def __init__(self, A=0.0, T3=0.0, T1=0.0, T2=None):
        """
        Parameters: *A*, *T3*, *T1*, *T2*. These must be entered as pure
        numbers with no attached dimensions.
        """
        self.A = A
        self.T3 = T3
        self.T1 = T1
        self.T2 = T2

    def get_yaml(self, out):
        troe = FlowMap([('A', self.A), ('T3', self.T3), ('T1', self.T1)])
        if self.T2 is not None:
            troe['T2'] = self.T2
        out['Troe'] = troe


class SRI(object):
    """ The SRI falloff function."""
    def __init__(self, A=0.0, B=0.0, C=0.0, D=None, E=None):
        """
        Parameters: *A*, *B*, *C*, *D*, *E*. These must be entered as
        pure numbers without attached dimensions.
        """
        self.A = A
        self.B = B
        self.C = C
        self.D = D
        self.E = E

    def get_yaml(self, out):
        sri = FlowMap([('A', self.A), ('B', self.B), ('C', self.C)])
        if self.D is not None and self.E is not None:
            sri['D'] = self.D
            sri['E'] = self.E
        out['SRI'] = sri


class Lindemann(object):
    """The Lindemann falloff function."""
    def get_yaml(self, out):
        pass


def convert(filename=None, outName=None, text=None):
    # Reset global state, in case cti2yaml is being used as a module and convert
    # is being called multiple times.
    units('m', 'kmol', 'kg', 's', 'J/kmol', 'J', 'Pa')
    standard_pressure(OneAtm)
    global _motz_wise
    _motz_wise = None
    _elements.clear()
    _species.clear()
    _speciesnames.clear()
    _phases.clear()
    _reactions.clear()
    _reactions['reactions'] = []

    if filename is not None:
        filename = os.path.expanduser(filename)
        base = os.path.basename(filename)
        root, _ = os.path.splitext(base)
        dataset(root)
    if outName is None and _name != 'noname':
        outName = _name + '.yaml'

    open_kw = {'encoding': 'latin-1'} if sys.version_info.major == 3 else {}
    open_mode = 'r' if sys.version_info.major == 3 else 'rU'
    try:
        if filename is not None:
            with open(filename, open_mode, **open_kw) as f:
                code = compile(f.read(), filename, 'exec')
        else:
            code = compile(text, '<string>', 'exec')
        exec(code)
    except SyntaxError as err:
        # Show more context than the default SyntaxError message
        # to help see problems in multi-line statements
        if filename:
            text = open(filename, open_mode).readlines()
        else:
            text = text.split('\n')
        _printerr('%s in "%s" on line %i:\n' % (err.__class__.__name__,
                                                err.filename,
                                                err.lineno))
        _printerr('|  Line |')
        for i in range(max(err.lineno-6, 0),
                       min(err.lineno+3, len(text))):
            _printerr('| % 5i |' % (i+1), text[i].rstrip())
            if i == err.lineno-1:
                _printerr(' '* (err.offset+9) + '^')
        _printerr()
        sys.exit(3)
    except Exception as err:
        import traceback

        if filename:
            text = open(filename, open_mode).readlines()
        else:
            text = text.split('\n')
            filename = '<string>'
        tb = traceback.extract_tb(sys.exc_info()[2])
        lineno = tb[-1][1]
        if tb[-1][0] == filename:
            # Error in input file
            _printerr('%s on line %i of %s:' % (err.__class__.__name__, lineno, filename))
            _printerr(err)
            _printerr('\n| Line |')

            for i in range(max(lineno-6, 0),
                           min(lineno+3, len(text))):
                if i == lineno-1:
                    _printerr('> % 4i >' % (i+1), text[i].rstrip())
                else:
                    _printerr('| % 4i |' % (i+1), text[i].rstrip())
        else:
            # Error in ctml_writer or elsewhere
            traceback.print_exc()

        sys.exit(4)

    # write the YAML file
    emitter = yaml.YAML()
    emitter.width = 70

    for name, cls in globals().items():
        if hasattr(cls, 'to_yaml'):
            emitter.register_class(cls)

    with open(outName, 'w') as dest:
        outputStarted = False

        # todo: header comment
        # todo: motz-wise

        outUnits = FlowMap([])
        if _umass != 'kg':
            outUnits['mass'] = _umass
        if _ulen != 'm':
            outUnits['length'] = _ulen
        if _utime != 's':
            outUnits['time'] = _utime
        if _upres != 'Pa':
            outUnits['pressure'] = _upres
        if _uenergy != 'J':
            outUnits['energy'] = _uenergy
        if _umol != 'kmol':
            outUnits['quantity'] = _umol
        if _ue != 'J/kmol':
            outUnits['activation-energy'] = _ue

        if outUnits:
            unitsMap = BlockMap([('units', outUnits)])
            emitter.dump(unitsMap, dest)
            outputStarted = True

        if _elements:
            elementsMap = BlockMap([('elements', _elements)])
            if outputStarted:
                elementsMap.yaml_set_comment_before_after_key('elements', before='\n')
            outputStarted = True
            emitter.dump(elementsMap, dest)

        if _phases:
            phasesMap = BlockMap([('phases', _phases)])
            if outputStarted:
                phasesMap.yaml_set_comment_before_after_key('phases', before='\n')
            outputStarted = True
            emitter.dump(phasesMap, dest)

        if _species:
            speciesMap = BlockMap([('species', _species)])
            if outputStarted:
                speciesMap.yaml_set_comment_before_after_key('species', before='\n')
            outputStarted = True
            emitter.dump(speciesMap, dest)

        for name, reactions in _reactions.items():
            if reactions:
                reactionsMap = BlockMap([(name, reactions)])
                if outputStarted:
                    reactionsMap.yaml_set_comment_before_after_key(name, before='\n')
                outputStarted = True
                emitter.dump(reactionsMap, dest)


def main():
    if len(sys.argv) == 1  or sys.argv[1] in ('-h', '--help'):
        print(__doc__)
        sys.exit(0)
    if len(sys.argv) not in (2,3):
        raise ValueError('Incorrect number of command line arguments.')
    convert(*sys.argv[1:])

if __name__ == "__main__":
    main()
