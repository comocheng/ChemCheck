import yaml
import math
import numpy as np
import matplotlib.pyplot as plt
import os
import linecache

def cp_calculate(T, T_mid, Nasa_poly_low, Nasa_poly_high):
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
    def check_negative_A_factor(self):
        arrhenius_reactions = []
        with open(self.path, 'r') as f:
            chem_data = yaml.load(f, Loader=yaml.FullLoader)
        reactions = chem_data['reactions']
        for r in reactions:
            if ('type', 'pressure-dependent-Arrhenius') in r.items():
                arrhenius_reactions.append(r)
        return arrhenius_reactions       

             
def err_line_without_comment(path, line_num):
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
    big_list.append(same_p_list)
    sum_of_selected_constants += len(same_p_list)
    rest_of_constants  = len(rate_constants) - sum_of_selected_constants
    if rest_of_constants > 0:
        same_p_list = list_of_rate_constants_with_same_pressure(rate_constants, sum_of_selected_constants, [])
        return bigger_list(rate_constants, same_p_list, big_list, sum_of_selected_constants)
    else:
        return big_list

class CheckNegativeA:
    def __init__(self, path):
        self.path = path
        
    def new_arrhenius_dict(self):
            a = {}
            b= {}
            #c = {}
            arrhenius_reactions = []
            #error_reactions = {}
            with open(self.path, 'r') as f:
                chem_data = yaml.load(f, Loader=yaml.FullLoader)
            reactions = chem_data['reactions']
            for r in reactions:
                if ('type', 'pressure-dependent-Arrhenius') in r.items():
                    arrhenius_reactions.append(r)
            #return arrhenius_reactions
            for reaction in arrhenius_reactions:
                rate_constants = reaction['rate-constants']
                reaction_equation = reaction['equation']
                b[reaction_equation] = rate_constants
            for i, k in b.items():
                same_p_list = list_of_rate_constants_with_same_pressure(k, 0, [])
                new_list_of_reaction_constants = bigger_list(k, same_p_list, [], 0)
                a[i] = new_list_of_reaction_constants
            return a

    def check_negative_A_factor(self, p):
        error_reactions = {}
        for equation, value in p.items():
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

    def check_sum_of_k(self, p, t):
        dict = {}
        for equation, value in p.items():
            dict[equation] = []
            for parameter_list in value:
                if len(parameter_list) != 1:
                    b = []
                    for rate_parameter in parameter_list:
                        rate_constant = rate_parameter['A'] * (t ** rate_parameter['b']) * math.exp(rate_parameter['Ea'] / (8314.46261815324 * t))
                        b.append(rate_constant)
                    c = sum(b)
                    if c <= 0:
                        dict[equation].append(parameter_list)
        error_equation_list = {}
        
        for key, value in dict.items():
            if dict[key] != []:
                error_equation_list[key] = value
        return error_equation_list
