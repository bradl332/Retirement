import numpy as np
import scipy
from scipy import optimize

from OLG_solve_steady import OLG_solve_steady
from funcs import prices, prodfunc

def calibration_util(x, global_dict):

    # x = beta, alpha

    print (x)

    if x[0] >= 0.995 or x[0] < 0.9:
        return 99999999

    if x[1] < 1.0:
        return 9999999

    args = {}
    args = dict.copy(global_dict)

    args['util_par']['beta'] = x[0]
    args['util_par']['alpha'] = x[1]

    K_Y_tgt = 3.0 / 5.0
    L_tgt = 0.3

    HH_res, K, L_dict = OLG_solve_steady(args)

    K_Y = K / prodfunc(K, L_dict, args)

    L = 0

    for tipo in L_dict.keys():
        L += args['pop'][tipo] * np.dot(args['pop_age'][tipo],HH_res[tipo][0])

    w,r = prices(K, L_dict, args)

    MOMKY = np.abs((K_Y - K_Y_tgt) / K_Y_tgt)
    MOML = np.abs((L - L_tgt) / L_tgt)

    print ('MOMKY',MOMKY)
    print ('MOML',MOML)

    return np.array([MOMKY, MOML])

# def calibration_eps(x, global_dict):
#
#     # x = beta, alpha
#
#     print (x)
#
#     args = {}
#     args = dict.copy(global_dict)
#
#     args['util_par']['beta'] = x[0]
#     args['util_par']['alpha'] = x[1]
#     args['eff']['RH'] = x[2]
#     args['eff']['AH'] = x[3]
#     args['eff']['PC'] = x[4]
#     args['eff']['RC'] = x[5]
#     args['eff']['AC'] = x[6]
#
#     K_Y_tgt = 3.0 / 5.0
#     L_tgt = 0.3
#
#     income_target_RH = 1.5
#     income_target_AH = 1.4
#     income_target_PC = 1.6
#     income_target_RC = 2.7
#     income_target_AC = 2.4
#
#     HH_res, K, L_dict = OLG_solve_steady(args)
#
#     Y = prodfunc(K, L_dict, args)
#
#     K_Y = K / Y
#
#     L = 0
#
#     for tipo in L_dict.keys():
#         L += args['pop'][tipo] * np.dot(args['pop_age'][tipo],HH_res[tipo][0])
#
#     w,r = prices(K, L_dict, args)
#
#     eff = args['eff']
#     pop_age = args['pop_age']
#
#     MOMKY = np.abs((K_Y - K_Y_tgt) / K_Y_tgt)
#     MOML = np.abs((L - L_tgt) / L_tgt)
#     MOME_RH = (((w['RH'] * eff['RH']) / w['PH']) - income_target_RH) / income_target_RH
#     MOME_AH = (((w['AH'] * eff['AH']) / w['PH']) - income_target_AH) / income_target_AH
#     MOME_PC = (((w['PC'] * eff['PC']) / w['PH']) - income_target_PC) / income_target_PC
#     MOME_RC = (((w['RC'] * eff['RC']) / w['PH']) - income_target_RC) / income_target_RC
#     MOME_AC = (((w['AC'] * eff['AC']) / w['PH']) - income_target_AC) / income_target_AC
#
#     MOME_AH = np.abs(MOME_AH)
#     MOME_RH = np.abs(MOME_RH)
#     MOME_PC = np.abs(MOME_PC)
#     MOME_RC = np.abs(MOME_RC)
#     MOME_AC = np.abs(MOME_AC)
#
#     print ('MOMKY',MOMKY)
#     print ('MOML',MOML)
#     print ('MOMAH',MOME_AH)
#     print ('MOMRH',MOME_RH)
#     print ('MOMRC',MOME_PC)
#     print ('MOMRC',MOME_RC)
#     print ('MOMAC',MOME_AC)
#
#     print ('result',MOMKY + MOML + MOME_AH + MOME_RH + MOME_PC + MOME_RC + MOME_AC)
#
#     return MOMKY + MOML + MOME_AH + MOME_RH + MOME_PC + MOME_RC + MOME_AC

def calibration_eps(x, global_dict):

    # x = beta, alpha

    print (x)

    if x[0] >= 0.995 or x[0] < 0.75:
        return 99999999

    if x[1] < 1.0:
        return 9999999

    if x[2] < 1.0:
        return 9999999

    if x[3] < 1.0:
        return 9999999

    if x[4] < 1.0:
        return 9999999

    if x[5] < 1.0:
        return 9999999

    if x[6] < 1.0:
        return 9999999

    if x[7] < 1.0:
        return 9999999

    if x[8] < 1.0:
        return 9999999

    args = {}
    args = dict.copy(global_dict)

    args['util_par']['beta'] = x[0]
    args['util_par']['alpha'] = x[1]
    args['eff']['RH'] = x[2]
    args['eff']['AH'] = x[3]
    args['eff']['PC'] = x[4]
    args['eff']['RC'] = x[5]
    args['eff']['ST'] = x[6]
    args['eff']['AI'] = x[7]
    args['eff']['AC'] = x[8]

    K_Y_tgt = 3.0 / 5.0
    L_tgt = 0.3

    income_target_RH = 1.74
    income_target_AH = 1.26
    income_target_PC = 1.70
    income_target_RC = 3.49
    income_target_ST = 2.91
    income_target_AI = 2.48
    income_target_AC = 2.25

    HH_res, K, L_dict = OLG_solve_steady(args)

    Y = prodfunc(K, L_dict, args)

    K_Y = K / Y

    L = 0

    for tipo in L_dict.keys():
        L += args['pop'][tipo] * np.dot(args['pop_age'][tipo],HH_res[tipo][0])

    w,r = prices(K, L_dict, args)

    eff = args['eff']
    pop_age = args['pop_age']

    MOMKY = np.abs((K_Y - K_Y_tgt) / K_Y_tgt)
    MOML = np.abs((L - L_tgt) / L_tgt)
    MOME_RH = (((w['RH'] * eff['RH']) / w['PH']) - income_target_RH) / income_target_RH
    MOME_AH = (((w['AH'] * eff['AH']) / w['PH']) - income_target_AH) / income_target_AH
    MOME_PC = (((w['PC'] * eff['PC']) / w['PH']) - income_target_PC) / income_target_PC
    MOME_RC = (((w['RC'] * eff['RC']) / w['PH']) - income_target_RC) / income_target_RC
    MOME_ST = (((w['ST'] * eff['ST']) / w['PH']) - income_target_ST) / income_target_ST
    MOME_AI = (((w['AI'] * eff['AI']) / w['PH']) - income_target_AI) / income_target_AI
    MOME_AC = (((w['AC'] * eff['AC']) / w['PH']) - income_target_AC) / income_target_AC

    MOME_AH = np.abs(MOME_AH)
    MOME_RH = np.abs(MOME_RH)
    MOME_PC = np.abs(MOME_PC)
    MOME_RC = np.abs(MOME_RC)
    MOME_ST = np.abs(MOME_ST)
    MOME_AI = np.abs(MOME_AI)
    MOME_AC = np.abs(MOME_AC)

    print ('MOMKY',MOMKY)
    print ('MOML',MOML)
    print ('MOMAH',MOME_AH)
    print ('MOMRH',MOME_RH)
    print ('MOMRC',MOME_PC)
    print ('MOMRC',MOME_RC)
    print ('MOMAC',MOME_AC)

    print ('result',MOMKY + MOML + MOME_AH + MOME_RH + MOME_PC + MOME_RC + MOME_AT + MOME_AI + MOME_AC)

    return MOMKY + MOML + MOME_AH + MOME_RH + MOME_PC + MOME_RC + MOME_AT + MOME_AI + MOME_AC
