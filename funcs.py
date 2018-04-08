import numpy as np
import scipy
from scipy import optimize

def prodfunc(K, L_types, global_dict):

    prod_par = global_dict['prod_par']

    L = agg_lab(L_types, global_dict)

    return prod_par['A'] * (K ** prod_par['theta']) * (L ** (1 - prod_par['theta']))

def agg_lab(L_types, global_dict):

    prod_par = global_dict['prod_par']
    lamda = prod_par['lamda']
    rho = prod_par['rho']

    L = 0

    # L_U = 0
    # L_S = 0
    #
    # for tipo in lamda['U'].keys():
    #     L_U += lamda['U'][tipo] * (L_types[tipo] ** ((rho['U'] - 1.0) / rho['U']))
    #
    # for tipo in lamda['S'].keys():
    #     L_S += lamda['S'][tipo] * (L_types[tipo] ** ((rho['S'] - 1.0) / rho['S']))
    #
    # L_U = L_U ** ((rho['U']) / (rho['U'] - 1.0))
    # L_S = L_S ** ((rho['S']) / (rho['S'] - 1.0))

    for tipo in L_types.keys():
        L += lamda['SU'][tipo] * (L_types[tipo] ** ((rho['SU'] - 1.0) / rho['SU']))

    # L = (lamda['SU']['U'] * (L_U ** ((rho['SU'] - 1.0) / rho['SU']))) + (lamda['SU']['S'] * (L_S ** ((rho['SU'] - 1.0) / rho['SU'])))

    L = L ** ((rho['SU']) / (rho['SU'] - 1.0))

    return L

def prices(K, L_types, global_dict):

    prod_par = global_dict['prod_par']

    wages = {}

    for t in L_types.keys():
        L_pe = dict.copy(L_types)
        L_ne = dict.copy(L_types)

        epsilon_L = L_pe[t] * 0.00001

        L_pe[t] = L_pe[t] + epsilon_L
        L_ne[t] = L_ne[t] - epsilon_L

        if L_pe[t] - epsilon_L == 0.0 or L_ne[t] + epsilon_L == 0.0:
            wages.update({t: 0})

        else:

            wages.update({t : (prodfunc(K, L_pe, global_dict) - prodfunc(K, L_ne, global_dict)) / (2 * epsilon_L)})

    epsilon_K = 0.00001 * K

    r = (prodfunc(K + epsilon_K, L_types, global_dict) - prodfunc(K - epsilon_K, L_types, global_dict)) / (2 * epsilon_K)

    r = r - prod_par['delta']

    return wages, r


def util(c, l, global_dict):

    util_par = global_dict['util_par']

    return np.log(c) + (util_par['alpha'] * np.log(1 - l))
