import upload.chemcheck as chemcheck
from django.test import TestCase
import yaml

class ChemicalTest(TestCase):
    def test_thermo_property_calculate(self):
        # This is the test to test the thermo properties' calculation results at mid temperature
        # the species has discontinuity in Cp, H, S
        species = {'name': 'C2H5CHO', 'composition': {'C': 3, 'H': 6, 'O': 1},
                   'thermo': {'model': 'NASA7', 'temperature-ranges': [273.15, 1500.0, 5000.0], 'data': [[7.6044596, -0.0086403564, 7.3930097e-05, -7.9687398e-08, 2.8004927e-11, -25489.789, -6.7643691], 
                   [3.3137982, 0.026619606, -1.0475596e-05, 1.8815334e-09, -1.276131e-13, -25459.603, 9.6608447]], 'note': 'T9/92'}, 
                   'transport': {'model': 'gas', 'geometry': 'nonlinear', 'well-depth': 357.0, 'diameter': 5.176, 'rotational-relaxation': 1.0, 'note': '=C4H8'}, 
                   'note': 'propionaldehyde'}
        T = species['thermo']['temperature-ranges']
        nasa_low = species['thermo']['data'][0]
        nasa_high = species['thermo']['data'][1]
        cp_1, cp_2 = chemcheck.cp_calculate(T[1], T[1], nasa_low, nasa_high)
        h_1, h_2 = chemcheck.h_calculate(T[1], T[1], nasa_low, nasa_high)
        s_1, s_2 = chemcheck.s_calculate(T[1], T[1], nasa_low, nasa_high)
        self.assertNotEqual(cp_1, cp_2)
        self.assertNotEqual(h_1, h_2)
        self.assertNotEqual(s_1, s_2)
