import numpy as np
import scipy
from scipy import optimize
from scipy.interpolate import make_interp_spline
import multiprocessing as mp
import sys
from scipy.interpolate import splrep, splev

from funcs import prodfunc, prices, util, agg_lab
from solv_func import agent_solve, agent_simul, para_agent, para_type_steady, gov_func

def OLG_solve_steady(global_dict):

    # Since this is partial equilbirium I can set prices to whatever the heck
    # I want. I use 1.0 and 0.5 because it makes plugging numbers into the algebra
    # a little easier

    eff = global_dict['eff']
    util_par =  global_dict['util_par']
    prod_par =  global_dict['prod_par']
    pop_age = global_dict['pop_age']
    pop = global_dict['pop']

    taxes = {}
    taxes.update({'labor tax': 0.0 * np.ones(global_dict['life'])})
    taxes.update({'capital tax': 0.0 * np.ones(global_dict['life'])})
    taxes.update({'transfers': 0.0 * np.ones(global_dict['life'])})

    ###################### Initial Guess of stuff ###############################

    per_length = 70.0 / global_dict['life']

    K_Y = 3.0 / per_length
    C_Y = 1 - (prod_par['delta'] * K_Y)

    C_Y_dict = {}
    L_dict = {}

    L_conso = (1.0 + ((util_par['alpha'] / (1.0 - prod_par['theta'])) * C_Y))
    L_conso = 1.0 / L_conso

    for tipo in global_dict['types_key']:

        L = np.sum(pop_age[tipo] * 0.3) * eff[tipo]
        L_dict[tipo] = L

    L = agg_lab(L_dict, global_dict)

    Y = (prod_par['A'] ** (1.0 / (1.0 - prod_par['theta']))) * L * (K_Y ** (prod_par['theta'] / (1.0 - prod_par['theta'])))

    K = K_Y * Y

    K_orig = np.copy(K)

    ###################### Starting the recursion #############################

    lam_K = 0.2
    iter_K = 0
    tol_K = 1 * (10 ** (- 4))
    close_K = False

    while close_K == False and iter_K < 100:

        w, r = prices(K, L_dict, global_dict)

        a_grid = np.linspace(0.0, global_dict['gen_par']['a grid scale'] * K, global_dict['points'])

        lam_L = 0.2
        iter_L = 0
        tol_L = 1.0 * (10 ** (-7))
        close_L = {}
        dist_L_dict = {}
        L_dict_new = {}

        asset_dict = {}

        for tipo in global_dict['types_key']:

            close_L.update({tipo: False})
            dist_L_dict.update({tipo: False})
            L_dict_new.update({tipo: False})

        while (any(list(close_L.values())) == False) and iter_L < 100:

            w, r = prices(K, L_dict, global_dict)

            L_dict_new = {}

            # Object args is global dict + stuff that can change on each recursion
            # version and V_nex are there to specify whether I want to solve the problem
            # of agents in retirement and also to provide a value function for death.

            object_args = dict.copy(global_dict)

            object_args.update({'a_grid': a_grid})
            object_args.update({'wages': w})
            object_args.update({'r': r})
            object_args.update({'V_nex': np.zeros(global_dict['points'])})
            object_args.update({'version': ''})
            object_args.update({'taxes': taxes})

            # For paralelization you can only pass one argument to the function
            # I pass a dictionary with the only thing changing is the type. For type
            # dependent parameters I can call them from within the solving code using
            # the type that is being solved for

            para_type_arg = []
            HH_results = {}

            for tipo in global_dict['types_key']:

                HH_results.update({tipo: None})
                asset_dict.update({tipo: None})

                type_object_args = dict.copy(object_args)
                type_object_args['types_key'] = tipo
                para_type_arg.append(type_object_args)

            pool = mp.Pool(len(global_dict['types_key']))

            results = list(pool.map(para_type_steady, para_type_arg))

            pool.close()
            pool.join()

            HH_res = {}

            for j in range(len(results)):
                tipo = results[j][-1]
                HH_results[tipo] = results[j][:-1]

            for tipo in global_dict['types_key']:

                ########### Need to check on labour supply #####################

                L_s = eff[tipo] * np.dot(pop_age[tipo], HH_results[tipo][0])

                L_dict_new[tipo] = L_s

                dist_L = np.abs((L_dict_new[tipo] - L_dict[tipo]) / L_dict[tipo])
                dist_L_dict[tipo] = dist_L

                close_L[tipo] = dist_L < tol_L

                L_dict[tipo] = (lam_L * L_dict_new[tipo]) + ((1.0 - lam_L) * L_dict[tipo])

                ########## Aggregating each types asset choices ################

                asset_tipo = np.dot(pop_age[tipo], HH_results[tipo][1]) - \
                            (pop_age[tipo][0] * HH_results[tipo][1][0]) + \
                            ((1.0 + global_dict['gen_par']['pop_growth']) * pop_age[tipo][0] * HH_results[tipo][1][0])

                asset_dict[tipo] = asset_tipo

            iter_L += 1

            # print ('L_dist',max(list(dist_L_dict.values())))

        K_new = sum(list(asset_dict.values()))

        dist_K = np.abs(K_new - K) / K

        print ('K dist',dist_K)

        close_K = dist_K < tol_K

        K = (lam_K * K_new) + ((1.0 - lam_K) * K)

        iter_K += 1

    print ('K dist',dist_K)

    return HH_results, K, L_dict
