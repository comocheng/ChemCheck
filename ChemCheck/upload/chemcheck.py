import yaml
import math
import numpy as np
import matplotlib.pyplot as plt
import os
import linecache
import cantera as ct

def cp_calculate(T, T_mid, Nasa_poly_low, Nasa_poly_high):
    """
    This function calculates NASA 7 polynomial Cp/R
    """
    if T < T_mid:
        cp_low = Nasa_poly_low[0] + Nasa_poly_low[1] * T + Nasa_poly_low[2] * T ** 2 + Nasa_poly_low[3] * T ** 3 + Nasa_poly_low[4] * T ** 4
        return cp_low
    elif T > T_mid:
        cp_high = Nasa_poly_high[0] + Nasa_poly_high[1] * T + Nasa_poly_high[2] * T ** 2 + Nasa_poly_high[3] * T ** 3 + Nasa_poly_high[4] * T ** 4
        return cp_high
    else:
        cp_low = Nasa_poly_low[0] + Nasa_poly_low[1] * T + Nasa_poly_low[2] * T ** 2 + Nasa_poly_low[3] * T ** 3 + Nasa_poly_low[4] * T ** 4
        cp_high = Nasa_poly_high[0] + Nasa_poly_high[1] * T + Nasa_poly_high[2] * T ** 2 + Nasa_poly_high[3] * T ** 3 + Nasa_poly_high[4] * T ** 4
        if cp_low - cp_high == 0:
            return cp_low
        else:
            return cp_low, cp_high

def h_calculate(T, T_mid, Nasa_poly_low, Nasa_poly_high):
    """
    This function calculates NASA 7 polynomial H/RT
    """
    if T < T_mid:
       h_low = Nasa_poly_low[0] + Nasa_poly_low[1] / 2 * T + Nasa_poly_low[2] / 3 * T ** 2 + Nasa_poly_low[3] / 4 * T ** 3 + Nasa_poly_low[4] / 5 * T ** 4 + Nasa_poly_low[5] / T
       return h_low
    elif T > T_mid:
        h_high = Nasa_poly_high[0] + Nasa_poly_high[1] / 2 * T + Nasa_poly_high[2] / 3 * T ** 2 + Nasa_poly_high[3] / 4 * T ** 3 + Nasa_poly_high[4] / 5 * T ** 4 + Nasa_poly_high[5] / T
        return h_high
    else:
        h_low = Nasa_poly_low[0] + Nasa_poly_low[1] / 2 * T + Nasa_poly_low[2] / 3 * T ** 2 + Nasa_poly_low[3] / 4 * T ** 3 + Nasa_poly_low[4] / 5 * T ** 4 + Nasa_poly_low[5] / T
        h_high = Nasa_poly_high[0] + Nasa_poly_high[1] / 2 * T + Nasa_poly_high[2] / 3 * T ** 2 + Nasa_poly_high[3] / 4 * T ** 3 + Nasa_poly_high[4] / 5 * T ** 4 + Nasa_poly_high[5] / T
        if h_low - h_high == 0:
            return h_low
        else:
            return h_low, h_high

def s_calculate(T, T_mid, Nasa_poly_low, Nasa_poly_high):
    """
    This function calculates NASA 7 polynomial S/R
    """
    if T < T_mid:
       s_low = Nasa_poly_low[0] * math.log(T) + Nasa_poly_low[1] * T + Nasa_poly_low[2] / 2 * T ** 2 + Nasa_poly_low[3] / 3 * T ** 3 + Nasa_poly_low[4] / 4 * T ** 4 + Nasa_poly_low[6]
       return s_low
    elif T > T_mid:
        s_high = Nasa_poly_high[0] * math.log(T) + Nasa_poly_high[1] * T + Nasa_poly_high[2] / 2 * T ** 2 + Nasa_poly_high[3] / 3 * T ** 3 + Nasa_poly_high[4] / 4 * T ** 4 + Nasa_poly_high[6]
        return s_high
    else:
        s_low = Nasa_poly_low[0] * math.log(T) + Nasa_poly_low[1] * T + Nasa_poly_low[2] / 2 * T ** 2 + Nasa_poly_low[3] / 3 * T ** 3 + Nasa_poly_low[4] / 4 * T ** 4 + Nasa_poly_low[6]
        s_high = Nasa_poly_high[0] * math.log(T) + Nasa_poly_high[1] * T + Nasa_poly_high[2] / 2 * T ** 2 + Nasa_poly_high[3] / 3 * T ** 3 + Nasa_poly_high[4] / 4 * T ** 4 + Nasa_poly_high[6]
        if s_low - s_high == 0:
            return s_low
        else:
            return s_low, s_high
            
def get_kf(rate_parameter, T):
    '''
    Calculate kf in m^3 / mol.s from Ea in cal/mol
    '''
    kf = rate_parameter['A'] * (T ** rate_parameter['b']) * math.exp(-rate_parameter['Ea'] / (1.985877534 * T)) / (1e6)
    return kf

class ChemicalError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)   

class ChemError:
    """
    This class checks the NASA Polynomial discontinuity and plot the figures
    for discontinuous species by matplotlib
    """
    def __init__(self, path, name):
        self.path = path
        self.name = name
        
    
    def check_continuity(self):
        discontinuous_species = []
        with open(self.path, 'r') as f:
            chem_data = yaml.load(f, Loader=yaml.FullLoader)
        for species in chem_data['species']:
            T_low = species['thermo']['temperature-ranges'][0]
            T_mid = species['thermo']['temperature-ranges'][1]
            #T_high = species['thermo']['temperature-ranges'][2]
            Nasa_poly_high = species['thermo']['data'][1]
            Nasa_poly_low = species['thermo']['data'][0]
            cp_low = Nasa_poly_low[0] + Nasa_poly_low[1] * T_mid + Nasa_poly_low[2] * T_mid ** 2 + Nasa_poly_low[3] * T_mid ** 3 + Nasa_poly_low[4] * T_mid ** 4
            cp_high = Nasa_poly_high[0] + Nasa_poly_high[1] * T_mid + Nasa_poly_high[2] * T_mid ** 2 + Nasa_poly_high[3] * T_mid ** 3 + Nasa_poly_high[4] * T_mid ** 4
            h_low = Nasa_poly_low[0] + Nasa_poly_low[1] / 2 * T_mid + Nasa_poly_low[2] / 3 * T_mid ** 2 + Nasa_poly_low[3] / 4 * T_mid ** 3 + Nasa_poly_low[4] / 5 * T_mid ** 4 + Nasa_poly_low[5] / T_mid
            h_high = Nasa_poly_high[0] + Nasa_poly_high[1] / 2 * T_mid + Nasa_poly_high[2] / 3 * T_mid ** 2 + Nasa_poly_high[3] / 4 * T_mid ** 3 + Nasa_poly_high[4] / 5 * T_mid ** 4 + Nasa_poly_high[5] / T_mid
            s_low = Nasa_poly_low[0] * math.log(T_mid) + Nasa_poly_low[1] * T_mid + Nasa_poly_low[2] / 2 * T_mid ** 2 + Nasa_poly_low[3] / 3 * T_mid ** 3 + Nasa_poly_low[4] / 4 * T_mid ** 4 + Nasa_poly_low[6]
            s_high = Nasa_poly_high[0] * math.log(T_mid) + Nasa_poly_high[1] * T_mid + Nasa_poly_high[2] / 2 * T_mid ** 2 + Nasa_poly_high[3] / 3 * T_mid ** 3 + Nasa_poly_high[4] / 4 * T_mid ** 4 + Nasa_poly_high[6]
            if abs(cp_high - cp_low) <= 0.1 and abs(h_low - h_high) <= 0.1 and abs(s_low - s_high) <= 0.1:
                pass
            else:
                image = '{}.png'.format(species['name'])
                img_path = os.path.join(os.path.split(self.path)[0], image)
                img_url = img_path.split('ChemCheck')[2]
                discontinuous_species.append(img_url)
                #gas = ct.Solution(self.path)
                #sp = gas.species(species['name'])
                T = np.linspace(T_low, 5000, 200)
                fig,ax = plt.subplots(1,3,figsize=(8,3.5))
                cp = [cp_calculate(tt, T_mid, Nasa_poly_low, Nasa_poly_high) for tt in T]
                h = [h_calculate(tt, T_mid, Nasa_poly_low, Nasa_poly_high) for tt in T]
                s = [s_calculate(tt, T_mid, Nasa_poly_low, Nasa_poly_high) for tt in T]
                print(species)
                ax[0].plot(T,cp)
                ax[0].set_title('$c_p/R$')
                ax[1].plot(T,h)
                ax[1].set_title('$h/RT$')
                ax[2].plot(T,s)
                ax[2].set_title('$s/R$')
                fig.tight_layout()
                fig.suptitle('{} thermo property'.format(species['name']), y=1)
                fig.savefig(img_path)
        return discontinuous_species


                # T_mid = symbols('x')
                # cp_low = Nasa_poly_low[0] + Nasa_poly_low[1] * T_mid + Nasa_poly_low[2] * T_mid ** 2 + Nasa_poly_low[3] * T_mid ** 3 + Nasa_poly_low[4] * T_mid ** 4
                # cp_high = Nasa_poly_high[0] + Nasa_poly_high[1] * T_mid + Nasa_poly_high[2] * T_mid ** 2 + Nasa_poly_high[3] * T_mid ** 3 + Nasa_poly_high[4] * T_mid ** 4
                # h_low = Nasa_poly_low[0] + Nasa_poly_low[1] / 2 * T_mid + Nasa_poly_low[2] / 3 * T_mid ** 2 + Nasa_poly_low[3] / 4 * T_mid ** 3 + Nasa_poly_low[4] / 5 * T_mid ** 4 + Nasa_poly_low[5] / T_mid
                # h_high = Nasa_poly_high[0] + Nasa_poly_high[1] / 2 * T_mid + Nasa_poly_high[2] / 3 * T_mid ** 2 + Nasa_poly_high[3] / 4 * T_mid ** 3 + Nasa_poly_high[4] / 5 * T_mid ** 4 + Nasa_poly_high[5] / T_mid
                # s_low = Nasa_poly_low[0] * ln(T_mid) + Nasa_poly_low[1] * T_mid + Nasa_poly_low[2] / 2 * T_mid ** 2 + Nasa_poly_low[3] / 3 * T_mid ** 3 + Nasa_poly_low[4] / 4 * T_mid ** 4 + Nasa_poly_low[6]
                # s_high = Nasa_poly_high[0] * ln(T_mid) + Nasa_poly_high[1] * T_mid + Nasa_poly_high[2] / 2 * T_mid ** 2 + Nasa_poly_high[3] / 3 * T_mid ** 3 + Nasa_poly_high[4] / 4 * T_mid ** 4 + Nasa_poly_high[6]
                # T_cp = solveset(Abs(cp_high - cp_low) < 1, T_mid, domain=S.Reals)
                # T_h = solveset(Abs(h_high - h_low) < 1, T_mid, domain=S.Reals)
                # T_s = solveset(Abs(s_high - s_low) < 1, T_mid, domain=S.Reals)
                # T_range = Interval(T_low, T_high)
                # T_continuous = T_cp.intersect(T_h).intersect(T_s).intersect(T_range)
                # T_continuouse = str(T_continuous)[str(T_continuous).find('Interval'):len(str(T_continuous))]
                # if not T_continuouse:
                #     print(species['name'])
                # else:
                #     print(species['name'], T_continuouse)

def err_line_without_comment(path, line_num):
    """
    This function helps locate the real postion where the thermo data of a
    species starts
    """
    err_line = linecache.getline(path, line_num)
    if len(err_line) != 0:
        if err_line[0] == '!' or err_line[-2] == '!':
            line_num += 1
            err_line = linecache.getline(path, line_num)
            return err_line_without_comment(path, line_num)
        else:
            return err_line, line_num
    else:
        return err_line

def list_of_rate_constants_with_same_pressure(rate_constants, index, same_pressure_list):
    """
    This helper function puts the arrhenius parameter under same pressure into one list
    """
    if index < len(rate_constants) - 1:
        if float(rate_constants[index]['P'].split()[0]) == float(rate_constants[index + 1]['P'].split()[0]):
            same_pressure_list.append(rate_constants[index])
            index += 1
            max_index = len(rate_constants) - 1
            if index == max_index:            
                same_pressure_list.append(rate_constants[index])
                return same_pressure_list
            else:
                return list_of_rate_constants_with_same_pressure(rate_constants, index, same_pressure_list)
        else:
            same_pressure_list.append(rate_constants[index])
            return same_pressure_list
    else:
        same_pressure_list.append(rate_constants[index])
        return same_pressure_list

def bigger_list(rate_constants, same_p_list, big_list, sum_of_selected_constants):
    """
    This helper function rearrange the thermo data for a species to make a bigger list 
    with sublists which contains arrhenius parameters of that species under the same pressure
    """
    big_list.append(same_p_list)
    sum_of_selected_constants += len(same_p_list)
    rest_of_constants  = len(rate_constants) - sum_of_selected_constants
    if rest_of_constants > 0:
        same_p_list = list_of_rate_constants_with_same_pressure(rate_constants, sum_of_selected_constants, [])
        return bigger_list(rate_constants, same_p_list, big_list, sum_of_selected_constants)
    else:
        return big_list

class CheckNegativeA:
    """
    This class check for negative A, negative sum of k for pressure dependent 
    reactions under same pressure, and negative sum of k for duplicate reactions 
    """
    def __init__(self, path):
        self.path = path
        with open(self.path, 'r') as f:
            self.chem_data = yaml.load(f, Loader=yaml.FullLoader)
        
    def new_arrhenius_dict(self):
        """
        This function rearrange the pressure dependent reactions 
        """
        new_arrhenius_dict = {}
        cons_list_of_1_rxn = {}
        arrhenius_reactions = []
        reactions = self.chem_data['reactions']
        for r in reactions:
            if ('type', 'pressure-dependent-Arrhenius') in r.items() and 'duplicate' not in r.keys():
                arrhenius_reactions.append(r)
        for reaction in arrhenius_reactions:
            rate_constants = reaction['rate-constants']
            reaction_equation = reaction['equation']
            cons_list_of_1_rxn[reaction_equation] = rate_constants
        for i, k in cons_list_of_1_rxn.items():
            same_p_list = list_of_rate_constants_with_same_pressure(k, 0, [])
            new_list_of_reaction_constants = bigger_list(k, same_p_list, [], 0)
            new_arrhenius_dict[i] = new_list_of_reaction_constants
        return new_arrhenius_dict

    def check_negative_A_factor(self, new_arrhenius_dict):
        """
        This function checks the negative A for pressure with only one set of arrhenius parameters
        """
        error_reactions = {}
        for equation, value in new_arrhenius_dict.items():
            for parameter_list in value:
                if len(parameter_list) == 1:
                    if parameter_list[0]['A'] <= 0:
                        if equation not in error_reactions.keys():
                            error_reactions[equation] = parameter_list
                        else:
                            error_reactions[equation].append(parameter_list)
                    else:
                        pass
        return error_reactions

    def check_sum_of_k(self, new_arrhenius_dict, t):
        """
        This fucntion checks the sum of k for pressure dependent reactions under 
        the pressure which has more than one set of arrhenius parameters
        """
        error_rxn_dict = {}
        for equation, value in new_arrhenius_dict.items():
            error_rxn_dict[equation] = []
            for parameter_list in value:
                if len(parameter_list) != 1:
                    k_list = []
                    for rate_parameter in parameter_list:
                        rate_constant = rate_parameter['A'] * (t ** rate_parameter['b']) * math.exp((-rate_parameter['Ea']) / (1.985877534 * t))
                        k_list.append(rate_constant)
                    sum_of_k = sum(k_list)
                    if sum_of_k <= 0:
                        #print(equation, sum_of_k)
                        error_rxn_dict[equation].append(parameter_list)
        error_equation_list = {}
        
        for key, value in error_rxn_dict.items():
            if error_rxn_dict[key] != []:
                error_equation_list[key] = value
        return error_equation_list
    
    def duplicate_reactions(self):
        """
        This function rearrange the pressure independent duplicate reactions to 
        input format for check_sum_of_k
        """
        duplicate_reactions = {}
        reactions = self.chem_data['reactions']
        for r in reactions:
            if ('duplicate', True) in r.items() and 'rate-constant' in r.keys():
                reaction_equation = r['equation']
                rate_constant = r['rate-constant']
                if reaction_equation not in duplicate_reactions.keys():
                    duplicate_reactions[reaction_equation] = [rate_constant]
                else:
                    duplicate_reactions[reaction_equation].append(rate_constant)
        for reaction_equation, rate_constant in duplicate_reactions.items():
            duplicate_reactions[reaction_equation] = [rate_constant]
        return duplicate_reactions
    
    def duplicate_reactions_multi_P(self):
        """
        This function rearrange the pressure dependent duplicate reactions to the input
        format for check_sum_of_k function
        """
        duplicate_reactions = {}
        reactions = self.chem_data['reactions']
        for r in reactions:
            if ('duplicate', True) in r.items() and 'rate-constants' in r.keys():
                reaction_equation = r['equation']
                rate_constants = r['rate-constants']
                if reaction_equation not in duplicate_reactions.keys():
                    duplicate_reactions[reaction_equation] = [rate_constants]
                else:
                    duplicate_reactions[reaction_equation].append(rate_constants)
        arrhenius_parameter_dict = {}
        new_arrhenius_parameter_dict = {}
        for rxn, rate_constants in duplicate_reactions.items():
            if len(rate_constants) == 1:
                for rate_constant in rate_constants:
                    same_p_list = list_of_rate_constants_with_same_pressure(rate_constant, 0, [])
                    new_list_of_reaction_constants = bigger_list(rate_constant, same_p_list, [], 0)
                    new_arrhenius_parameter_dict[rxn] = new_list_of_reaction_constants
            else:
                arrhenius_parameter_dict[rxn] = []
                for list in rate_constants:
                    for parameter_dict in list:
                        arrhenius_parameter_dict[rxn].append(parameter_dict)
        for parameter_list in arrhenius_parameter_dict.values():
            for idx_of_para in range(len(parameter_list)):
                min_idx = idx_of_para
                for idx_of_rest_para in range(idx_of_para + 1, len(parameter_list)):
                    if parameter_list[min_idx]['P'].split()[0] > parameter_list[idx_of_rest_para]['P'].split()[0]:
                        min_idx = idx_of_rest_para
                parameter_list[idx_of_para], parameter_list[min_idx] = parameter_list[min_idx], parameter_list[idx_of_para]
        for i, k in arrhenius_parameter_dict.items():
                    same_p_list = list_of_rate_constants_with_same_pressure(k, 0, [])
                    new_list_of_reaction_constants = bigger_list(k, same_p_list, [], 0)
                    new_arrhenius_parameter_dict[i] = new_list_of_reaction_constants
        return new_arrhenius_parameter_dict



class check_collision_violation(CheckNegativeA):
    def __init__(self, path):
        self.path = path
        self.chem_elements = {'H': 1.00811, 'He': 4.002602, 'Li':6.938, 'Be':9.0121831,'B':10.806, 
                              'C': 12.0096, 'N': 14.00728, 'O': 15.99903, 'F': 18.998403163, 'Ne': 20.1797, 
                              'Si': 28.084,'P': 30.973761998, 'S': 32.059, 'Cl': 35.446, 'Ar': 39.948, 
                              'Br': 79.904, 'I': 126.90447, 'Kr': 83.8, 'Xe': 131.29}

        with open(self.path, 'r') as f:
            self.chem_data = yaml.load(f, Loader=yaml.FullLoader)
        if 'elements' in self.chem_data.keys():
            isotope_list = self.chem_data['elements']
            for isotope in isotope_list:
                symbol = isotope['symbol']
                atomic_weight = isotope['atomic-weight']
                try:
                    self.chem_elements[symbol]
                    raise Exception('isotope name is same with the element, please change')
                except KeyError:
                    self.chem_elements[symbol] = atomic_weight
        
    def get_reduced_mass(self, equation, reverse=False):
        '''
        reduced mass unit is in kg
        '''
        if not reverse:
            reactants = self.get_reactants(equation)
            for i, reactant in enumerate(reactants):
                if reactant == 'NO':
                    reactants[i] = False
            species_list = self.chem_data['species']
            reactant1 = 0
            reactant2 = 0 
            for species in species_list:
                if species['name'] == reactants[0]:
                    composition = species['composition']
                    for element, num in composition.items():
                        reactant1 += self.chem_elements[element] * num * (1e-3)
                if species['name'] == reactants[1]:
                    composition = species['composition']
                    for element, num in composition.items():
                        reactant2 += self.chem_elements[element] * num * (1e-3)
            if reactant1 == 0:
                print(reactants, 'reactants', equation)
                raise ChemicalError('No species data for {}'.format(reactants[0]))
            elif reactant2 == 0:
                print(reactants, 'reactants', equation)
                raise ChemicalError('No species data for {}'.format(reactants[1]))
            reduced_mass = reactant1 * reactant2 / (reactant1 + reactant2)
            return reduced_mass
        else:
            products = self.get_products(equation)
            for i, product in enumerate(products):
                if product == 'NO':
                    products[i] = False
            species_list = self.chem_data['species']
            product1 = 0
            product2 = 0 
            for species in species_list:
                if species['name'] == products[0]:
                    composition = species['composition']
                    for element, num in composition.items():
                        product1 += self.chem_elements[element] * num * (1e-3)
                if species['name'] == products[1]:
                    composition = species['composition']
                    for element, num in composition.items():
                        product2 += self.chem_elements[element] * num * (1e-3)
            if product1 == 0:
                print(products, 'products', equation)
                raise ChemicalError('No species data for {}'.format(products[0]))
            elif product2 == 0:
                print(products, 'products', equation)
                raise ChemicalError('No species data for {}'.format(products[1]))
            reduced_mass = product1 * product2 / (product1 + product2)
            return reduced_mass
        
    def get_epsilon_sigma(self, equation, reverse=False):
        epsilon = []
        sigma = []
        if reverse:
            products = self.get_products(equation)
            for product in products:
                for spc in self.chem_data['species']:
                    if spc['name'] == product:
                        well_depth = spc['transport']['well-depth'] * 1.3807e-23 * 6.02214179e23 * 1000
                        diameter = spc['transport']['diameter'] * (1e-10)
                        epsilon.append(well_depth)
                        sigma.append(diameter)
            mean_epsilon = (epsilon[0] * epsilon[1]) ** (1 / len(products))
            mean_sigma = sum(sigma) / len(products)
            return mean_epsilon, mean_sigma
        else:
            reactants = self.get_reactants(equation)
            for reactant in reactants:
                for spc in self.chem_data['species']:
                    if spc['name'] == reactant:
                        well_depth = spc['transport']['well-depth'] * 1.3807e-23 * 6.02214179e23
                        diameter = spc['transport']['diameter'] * (1e-10)
                        epsilon.append(well_depth)
                        sigma.append(diameter)
            mean_epsilon = (epsilon[0] * epsilon[1]) ** (1 / len(reactants))
            mean_sigma = sum(sigma) / len(reactants)
            return mean_epsilon, mean_sigma

    def cal_collision_limit(self,T, equation, reverse=False):
        '''
        The collision limit unit is in m^3 /mol.s
        '''
        if not reverse:
            reduced_mass = self.get_reduced_mass(equation)
            epsilon, sigma = self.get_epsilon_sigma(equation)
            Tr = T * (1.3806504e-23) * (6.02214179e23) / epsilon
            reduced_coll_integral = 1.16145 * Tr ** (-0.14874) + 0.52487 * math.exp(-0.7732 * Tr) + 2.16178 * math.exp(
                                    -2.437887 * Tr)
            k_coll = (math.sqrt(8 * math.pi * (1.3806504e-23) * T * (6.02214179e23) / reduced_mass) * sigma ** 2
                  * reduced_coll_integral * (6.02214179e23))
            return k_coll
        else:
            reduced_mass = self.get_reduced_mass(equation, reverse=True)
            epsilon, sigma = self.get_epsilon_sigma(equation, reverse=True)
            Tr = T * (1.3806504e-23) * (6.02214179e23) / epsilon
            reduced_coll_integral = 1.16145 * Tr ** (-0.14874) + 0.52487 * math.exp(-0.7732 * Tr) + 2.16178 * math.exp(
                                    -2.437887 * Tr)
            k_coll = (math.sqrt(8 * math.pi * (1.3806504e-23) * T * (6.02214179e23) / reduced_mass) * sigma ** 2
                     * reduced_coll_integral * (6.02214179e23))
            return k_coll
    
    def get_reactants(self, equation):
        reactants = []
        stoi_num = 1
        eqn_eles = equation.split()
        for idx, eqn_ele in enumerate(eqn_eles):
            if eqn_ele == '<=>' or eqn_ele == '=>' or eqn_ele == '=':
                first_half = eqn_eles[0:idx]
                for idx, reactant in enumerate(first_half):
                    if reactant.isdigit() == True:
                        stoi_num = int(reactant)
                    if reactant.isdigit() == False and reactant != '+':
                        separate_reactants = [first_half[idx]] * stoi_num
                        reactants.extend(separate_reactants)
                        stoi_num = 1
        for i, reactant in enumerate(reactants):
            if reactant == 'NO':
                reactants[i] = False      
        return reactants
    
    def get_products(self, equation):
        products = []
        stoi_num = 1
        eqn_eles = equation.split()
        for idx, eqn_ele in enumerate(eqn_eles):
            if eqn_ele == '<=>' or eqn_ele == '=>' or eqn_ele == '=':
                second_half = eqn_eles[idx + 1:len(eqn_eles)]
                for idx, product in enumerate(second_half):
                    if product.isdigit() == True:
                        stoi_num = int(product)
                    if product.isdigit() == False and product != '+':
                        separate_products = [second_half[idx]] * stoi_num
                        products.extend(separate_products)
                        stoi_num = 1
                    # if product == '+':
                    #     del products[idx]
                    # elif product.isdigit() == True:
                    #     separate_products = [products[idx + 1]] * int(product)
                    #     del products[idx:idx + 2]
                    #     products.extend(separate_products)
        for i, product in enumerate(products):
            if product == 'NO':
                products[i] = False
        return products
    
    def get_equilibrium_constant(self, reactants, products, T):
        reactants_H = 0
        products_H = 0
        reactants_S = 0
        products_S = 0
        for reactant in reactants:
            for species in self.chem_data['species']:
                if species['name'] == reactant:
                    T_mid = species['thermo']['temperature-ranges'][1]
                    Nasa_poly_high = species['thermo']['data'][1]
                    Nasa_poly_low = species['thermo']['data'][0]
                    h = h_calculate(T, T_mid, Nasa_poly_low, Nasa_poly_high)
                    s = s_calculate(T, T_mid, Nasa_poly_low, Nasa_poly_high)
                    if type(h) is tuple:
                        reactants_H += h[1] 
                    else:
                        reactants_H += h 
                    if type(s) is tuple:
                        reactants_S += s[1] 
                    else:
                        reactants_S += s 
        for product in products:
            for species in self.chem_data['species']:
                if species['name'] == product:
                    T_mid = species['thermo']['temperature-ranges'][1]
                    Nasa_poly_high = species['thermo']['data'][1]
                    Nasa_poly_low = species['thermo']['data'][0]
                    h = h_calculate(T, T_mid, Nasa_poly_low, Nasa_poly_high)
                    s = s_calculate(T, T_mid, Nasa_poly_low, Nasa_poly_high)
                    if type(h) is tuple:
                        products_H += h[1]
                    else:
                        products_H += h 
                    if type(s) is tuple:
                        products_S += s[1] 
                    else:
                        products_S += s 
        Gibbs_free_energy = (products_H - reactants_H) - (products_S - reactants_S) #delta G / RT
        equilibrium_constant = math.exp(-Gibbs_free_energy)
        return equilibrium_constant    
    
    def get_elementary_rxns(self):
        elementary_rxns = []
        reactions = self.chem_data['reactions']
        for rxn in reactions:
            if list(rxn.keys()) == ['equation', 'rate-constant'] or list(rxn.keys()) == ['equation', 'rate-constant', 'note']:
                elementary_rxns.append(rxn)
        return elementary_rxns
                
        
    def check_collision_violation(self, T, P):
        try:        
            gas = ct.Solution(self.path)
            gas.TP = T, P
            rxns = self.chem_data['reactions']
            violation_check = {}
            for index, rxn in enumerate(rxns):
                eqn = rxn['equation']
                reactants = self.get_reactants(eqn)
                products = self.get_products(eqn)
                if len(reactants) == 2 and 'falloff' not in rxn.values() and 'three-body' not in rxn.values() and 'M' not in reactants and '(+M)' not in reactants:
                    collision_limit = self.cal_collision_limit(T, eqn)
                    kf = gas.forward_rate_constants[index] / (1e3)
                    if collision_limit < kf:
                        violation_factor = kf / collision_limit
                        if eqn in violation_check.keys():
                            violation_check[eqn]['forward violation factor'] = violation_factor
                            violation_check[eqn]['forward collision limit'] = str(collision_limit) + ' m3/mol.s'
                            violation_check[eqn]['forward rate coefficient'] = str(kf) + + ' m3/mol.s'
                        else:
                            violation_check[eqn] = {'forward violation factor': violation_factor, 'forward collision limit': str(collision_limit) + ' m3/mol.s',
                                                'forward rate coefficient': str(kf) + ' m3/mol.s'}
                if len(products) == 2 and 'falloff' not in rxn.values() and 'three-body' not in rxn.values() and 'M' not in products and '(+M)' not in products:
                    collision_limit = self.cal_collision_limit(T, eqn, reverse=True)
                    kr = gas.reverse_rate_constants[index] / (1e3)
                    if collision_limit < kr:
                        violation_factor = kr / collision_limit
                        if eqn in violation_check.keys():
                            violation_check[eqn]['reverse violation factor'] = violation_factor
                            violation_check[eqn]['reverse collision limit'] = str(collision_limit) + ' m3/mol.s'
                            violation_check[eqn]['reverse rate coeffcient'] = str(kr) + ' m3/mol.s'
                        else:    
                            violation_check[eqn] = {'reverse violation factor': violation_factor, 'reverse collision limit': str(collision_limit) + ' m3/mol.s',
                                                    'reverse rate coefficient': str(kr) + ' m3/mol.s'}
            return violation_check

        except ct.CanteraError:
            P = P / 1e5
            violation_check = {}
            elementary_rxns = self.get_elementary_rxns()
            duplicate_rxns = self.duplicate_reactions()
            duplicate_rxns_multi_P = self.duplicate_reactions_multi_P()
            pdep_rxns = self.new_arrhenius_dict()
            
            for rxn in elementary_rxns:
                eqn = rxn['equation']
                reactants = self.get_reactants(eqn)
                products = self.get_products(eqn)
                if len(reactants) == 2 and 'falloff' not in rxn.values() and 'three-body' not in rxn.values() and 'M' not in reactants and '(+M)' not in reactants:
                    collision_limit = self.cal_collision_limit(T, eqn)
                    rate_parameter = rxn['rate-constant']
                    kf = get_kf(rate_parameter, T)
                    if collision_limit < kf:
                        violation_factor = kf / collision_limit
                        if eqn in violation_check.keys():
                            violation_check[eqn]['forward violation factor'] = violation_factor
                            violation_check[eqn]['forward collision limit'] = str(collision_limit) + ' m3/mol.s'
                            violation_check[eqn]['forward rate coefficient'] = str(kf) + ' m3/mol.s'
                        else:
                            violation_check[eqn] = {'forward violation factor': violation_factor, 'forward collision limit': str(collision_limit) + ' m3/mol.s',
                                                'forward rate coefficient': str(kf) + ' m3/mol.s'}
                if len(products) == 2 and 'falloff' not in rxn.values() and 'three-body' not in rxn.values() and 'M' not in products and '(+M)' not in products:
                    collision_limit = self.cal_collision_limit(T, eqn, reverse=True)
                    equilibrium_constant = self.get_equilibrium_constant(reactants, products, T)
                    rate_parameter = rxn['rate-constant']
                    if '=>' not in eqn:
                        kf = get_kf(rate_parameter, T)
                    else:
                        kf = 0
                    kr = kf / equilibrium_constant
                    if collision_limit < kr:
                        violation_factor = kr / collision_limit
                        if eqn in violation_check.keys():
                            violation_check[eqn]['reverse violation factor'] = violation_factor
                            violation_check[eqn]['reverse collision limit'] = str(collision_limit) + ' m3/mol.s'
                            violation_check[eqn]['reverse rate coefficient'] = str(kr) + ' m3/mol.s'
                        else:
                            violation_check[eqn] = {'reverse violation factor': violation_factor, 'reverse collision limit': str(collision_limit) + ' m3/mol.s',
                                                'reverse rate coefficient': str(kr) + ' m3/mol.s'}
            
            for eqn, arr_parameters in duplicate_rxns.items():
                reactants = self.get_reactants(eqn)
                products = self.get_products(eqn)
                if len(reactants) == 2 and 'falloff' not in rxn.values() and 'three-body' not in rxn.values() and 'M' not in reactants and '(+M)' not in reactants:
                    collision_limit = self.cal_collision_limit(T, eqn)
                    kf_list = []
                    for arr_parameter in arr_parameters[0]:
                        kf = get_kf(arr_parameter, T)
                        kf_list.append(kf)
                    kf_dup = sum(kf_list)
                    if collision_limit < kf_dup:
                        violation_factor = kf_dup / collision_limit
                        violation_check[eqn] = {'type': 'duplicate reaction', 'forward violation factor': violation_factor, 'forward collision limit': str(collision_limit) + ' m3/mol.s',
                                                'forward rate coefficient': str(kf_dup) + ' m3/mol.s'}
                if len(products) == 2 and 'falloff' not in rxn.values() and 'three-body' not in rxn.values() and 'M' not in products and '(+M)' not in products:
                    collision_limit = self.cal_collision_limit(T, eqn, reverse=True)
                    equilibrium_constant = self.get_equilibrium_constant(reactants, products, T)
                    if '=>' not in eqn:
                        for arr_parameter in arr_parameters[0]:
                            kf = get_kf(arr_parameter, T)
                            kf_list.append(kf)
                        kf_dup = sum(kf_list)
                    else:
                        kf_dup = 0
                    kr_dup = kf_dup / equilibrium_constant
                    if collision_limit < kr_dup:
                        violation_factor = kr_dup / collision_limit
                        if eqn in violation_check.keys():
                            violation_check[eqn]['reverse violation factor'] = violation_factor
                            violation_check[eqn]['reverse collision limit'] = str(collision_limit) + ' m3/mol.s'
                            violation_check[eqn]['reverse rate coeffcient'] = str(kr_dup) + ' m3/mol.s'
                        else:    
                            violation_check[eqn] = {'type': 'duplicate reaction', 'reverse violation factor': violation_factor, 'reverse collision limit': str(collision_limit) + ' m3/mol.s',
                                                    'reverse rate coefficient': str(kr_dup) + ' m3/mol.s'} 
            
            for eqn, arr_parameters in pdep_rxns.items():
                reactants = self.get_reactants(eqn)
                products = self.get_products(eqn)
                if len(reactants) == 2 and 'falloff' not in rxn.values() and 'three-body' not in rxn.values() and 'M' not in reactants and '(+M)' not in reactants:
                    collision_limit = self.cal_collision_limit(T, eqn)
                    kf_pdep = self.cal_pdep_rate(arr_parameters, T, P)
                    if kf_pdep > collision_limit:
                        violation_factor = kf_pdep / collision_limit
                        violation_check[eqn] = {'type': 'pressure dependent reaction', 'forward violation factor': violation_factor, 'collision limit': str(collision_limit) + ' m3/mol.s',
                                                'forward rate coefficient': str(kf_pdep) + ' m3/mol.s'}
                if len(products) == 2 and 'falloff' not in rxn.values() and 'three-body' not in rxn.values() and 'M' not in products and '(+M)' not in products:
                    collision_limit = self.cal_collision_limit(T, eqn, reverse=True)
                    equilibrium_constant = self.get_equilibrium_constant(reactants, products, T)
                    if '=>' not in eqn:
                        kf_pdep = self.cal_pdep_rate(arr_parameters, T, P)
                    else:
                        kf_pdep = 0
                    kr_pdep = kf_pdep / equilibrium_constant
                    if kr_pdep > collision_limit:
                        violation_factor = kr_pdep / collision_limit
                        if eqn in violation_check.keys():
                            violation_check[eqn]['reverse violation factor'] = violation_factor
                            violation_check[eqn]['reverse collision limit'] = str(collision_limit) + ' m3/mol.s'
                            violation_check[eqn]['reverse rate coeffcient'] = str(kr_pdep) + ' m3/mol.s'
                        else:
                            violation_check[eqn] = {'type': 'pressure dependent reaction', 'reverse violation factor': violation_factor, 'collision limit': str(collision_limit) + ' m3/mol.s',
                                                    'reverse rate coefficient': str(kr_pdep) + ' m3/mol.s'}
            
            for eqn, arr_parameters in duplicate_rxns_multi_P.items():
                reactants = self.get_reactants(eqn)
                products = self.get_products(eqn)
                if len(reactants) == 2 and 'falloff' not in rxn.values() and 'three-body' not in rxn.values() and 'M' not in reactants and '(+M)' not in reactants:
                    collision_limit = self.cal_collision_limit(T, eqn)
                    kf_pdep = self.cal_pdep_rate(arr_parameters, T, P)
                    if kf_pdep > collision_limit:
                        violation_factor = kf_pdep / collision_limit
                        violation_check[eqn] = {'type': 'duplicate pressure dependent reaction', 'forward violation factor': violation_factor, 'collision limit': str(collision_limit) + ' m3/mol.s',
                                                'forward rate coefficient': str(kf_pdep) + ' m3/mol.s'}
                if len(products) == 2 and 'falloff' not in rxn.values() and 'three-body' not in rxn.values() and 'M' not in products and '(+M)' not in products:
                    collision_limit = self.cal_collision_limit(T, eqn, reverse=True)
                    equilibrium_constant = self.get_equilibrium_constant(reactants, products, T)
                    if '=>' not in eqn:
                        kf_pdep = self.cal_pdep_rate(arr_parameters, T, P)
                    else:
                        kf_pdep = 0
                    kr_pdep = kf_pdep / equilibrium_constant
                    if kr_pdep > collision_limit:
                        violation_factor = kr_pdep / collision_limit
                        if eqn in violation_check.keys():
                            violation_check[eqn]['reverse violation factor'] = violation_factor
                            violation_check[eqn]['reverse collision limit'] = str(collision_limit) + ' m3/mol.s'
                            violation_check[eqn]['reverse rate coeffcient'] = str(kr_pdep) + ' m3/mol.s'
                        else:
                            violation_check[eqn] = {'type': 'duplicate pressure dependent reaction', 'reverse violation factor': violation_factor, 'collision limit': str(collision_limit) + ' m3/mol.s',
                                                    'reverse rate coefficient': str(kr_pdep) + ' m3/mol.s'}
                    
            return violation_check
    
    def cal_pdep_rate(self, pdep_rxns:list, T, P)->float:
        small_p = []
        large_p = []
        for arr_parameters in pdep_rxns:
            p = float(arr_parameters[0]['P'].split()[0])
            if p == P and len(arr_parameters) == 1:
                rate_parameters = arr_parameters[0]
                kf = get_kf(rate_parameters, T)
                return kf
            elif p == P and len(arr_parameters) > 1:
                kf_list = []
                for arr_parameter in arr_parameters:
                    kf = get_kf(arr_parameter, T)
                    kf_list.append(kf)
                kf_sum = sum(kf_list)
                return kf_sum  
            elif p < P:
                small_p.append(p)
            elif p > P:
                large_p.append(p)
        if len(large_p) != 0 and len(small_p) != 0:
            p1 = small_p[-1]
            p2 = large_p[0]
            for i, arr_parameters in enumerate(pdep_rxns):
                p = float(arr_parameters[0]['P'].split()[0])
                if p == p1:
                    k1_arr = pdep_rxns[i]
                    k1_list = []
                    for arr_parameter in k1_arr:
                        kf = get_kf(arr_parameter, T)
                        k1_list.append(kf)
                    kf1 = sum(k1_list)            
                elif p == p2:
                    k2_arr = pdep_rxns[i]
                    k2_list = []
                    for arr_parameter in k2_arr:
                        kf = get_kf(arr_parameter, T)
                        k2_list.append(kf)
                    kf2 = sum(k2_list)
            try:
                log_k_P = math.log(kf1) + (math.log(kf2) - math.log(kf1)) * (math.log(P / p1) / math.log(p2 / p1))
                k_P = math.exp(log_k_P)
            except ValueError:
                pass
            return k_P
        elif len(large_p) == 0:
            parameters = arr_parameters[-1]
            if len(parameters) == 1:
                rate_parameters = parameters[0]
                kf = get_kf(rate_parameters, T)
                return kf
            if len(parameters) > 1:
                kf_list = []
                for arr_parameter in arr_parameters:
                    kf = get_kf(arr_parameter, T)
                    kf_list.append(kf)
                kf_sum = sum(kf_list)
                return kf_sum
        elif len(small_p) == 0:
            parameters = arr_parameters[0]
            if len(parameters) == 1:
                rate_parameters = parameters[0]
                kf = get_kf(rate_parameters, T)
                return kf
            if len(parameters) > 1:
                kf_list = []
                for arr_parameter in arr_parameters:
                    kf = get_kf(arr_parameter, T)
                    kf_list.append(kf)
                kf_sum = sum(kf_list)
                return kf_sum
