import yaml
import math
import numpy as np
import matplotlib.pyplot as plt
import os
import linecache

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
            #T_low = species['thermo']['temperature-ranges'][0]
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
                T = np.linspace(500, 3000, 200)
                fig,ax = plt.subplots(1,3,figsize=(8,3.5))
                cp = [cp_calculate(tt, T_mid, Nasa_poly_low, Nasa_poly_high) for tt in T]
                h = [h_calculate(tt, T_mid, Nasa_poly_low, Nasa_poly_high) for tt in T]
                s = [s_calculate(tt, T_mid, Nasa_poly_low, Nasa_poly_high) for tt in T]
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
            if ('type', 'pressure-dependent-Arrhenius') in r.items():
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
                        rate_constant = rate_parameter['A'] * (t ** rate_parameter['b']) * math.exp(rate_parameter['Ea'] / (1.985877534 * t))
                        k_list.append(rate_constant)
                    sum_of_k = sum(k_list)
                    if sum_of_k <= 0:
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
