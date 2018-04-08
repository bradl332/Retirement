import numpy as np
import scipy
from scipy import optimize
from scipy.interpolate import make_interp_spline
import multiprocessing as mp
import pickle

from OLG_solve_steady import OLG_solve_steady
from calibration import calibration_util, calibration_eps
from funcs import prodfunc, prices, util, agg_lab

a = np.linspace(0,10,9)

work = 8
retire = 0
life = work + retire
na = 1200

util_par = {}
util_par.update({'beta': 0.925})
util_par.update({'alpha': 1.779})

rho = {}
rho.update({'SU': 2.6})
rho.update({'U': 5.6})
rho.update({'S': 20.0})

lamda = {}
lamda_su = {}
lamda_u = {}
lamda_s = {}

lamda_su.update({'S': 0.5})
lamda_su.update({'U': 0.5})

lamda_u.update({'PH': 0.28})
lamda_u.update({'RH': 0.16})
lamda_u.update({'AH': 0.35})
lamda_u.update({'PC': 0.21})

lamda_s.update({'RC': 0.51})
lamda_s.update({'AC': 0.49})

lamda.update({'SU': lamda_su})
lamda.update({'U': lamda_u})
lamda.update({'S': lamda_s})

prod_par = {}
prod_par.update({'theta': 0.36})

per_length = 40.0 / life
delta = (prod_par['theta'] * (1.0 / (3.0 / per_length)) - (((1.05) ** (per_length)) - 1)) + 0.2
prod_par.update({'delta': delta})
prod_par.update({'A': 1.0})
prod_par.update({'lamda': lamda})
prod_par.update({'rho': rho})

gen_par = {}
gen_par.update({'bequest': 0.00})
gen_par.update({'a grid scale': 3.0})

eff = {}
eff.update({'PH': 1.0})
eff.update({'RH': 1.5})
eff.update({'AH': 1.4})
eff.update({'PC': 1.6})
eff.update({'RC': 2.7})
eff.update({'AC': 2.4})

pop_age = {}
pop_age.update({'PH': (1.0 / life) * np.ones(life)})
pop_age.update({'RH': (1.0 / life) * np.ones(life)})
pop_age.update({'AH': (1.0 / life) * np.ones(life)})
pop_age.update({'PC': (1.0 / life) * np.ones(life)})
pop_age.update({'RC': (1.0 / life) * np.ones(life)})
pop_age.update({'AC': (1.0 / life) * np.ones(life)})

pop = {}

pop.update({'PH': 0.10})
pop.update({'RH': 0.005})
pop.update({'AH': 0.505})
pop.update({'PC': 0.02})
pop.update({'RC': 0.05})
pop.update({'AC': 0.32})

types = len(list(eff.keys()))
types_key = list(eff.keys())

global_dict = {}
global_dict.update({'util_par': util_par})
global_dict.update({'prod_par': prod_par})
global_dict.update({'eff': eff})
global_dict.update({'pop_age': pop_age})
global_dict.update({'pop': pop})
global_dict.update({'life': life})
global_dict.update({'work': work})
global_dict.update({'retire': retire})
global_dict.update({'types_key': types_key})
global_dict.update({'points': na})
global_dict.update({'gen_par': gen_par})

HH_res_baseline, K_baseline, L_dict = OLG_solve_steady(global_dict)

w_baseline,r_baseline = prices(K_baseline, L_dict, global_dict)

# util_baseline = {}
# util_baseline['AH'] = np.sum(util(HH_res['AH'][2], HH_res['AH'][0], global_dict))
# util_baseline['AC'] = np.sum(util(HH_res['AC'][2], HH_res['AC'][0], global_dict))

L_baseline = dict.copy(L_dict)

Y_baseline = prodfunc(K_baseline,L_baseline, global_dict)
L_agg_baseline = agg_lab(L_baseline, global_dict)

print ('Baseline complete')
print ('HH_res',HH_res_baseline)
print ('a_max',K_baseline * gen_par['a grid scale'])
print ('#############################################')
print ('RC Balance?')
C = 0

for tipo in global_dict['types_key']:
    C += pop[tipo] * np.dot(pop_age[tipo], HH_res_baseline[tipo][2])

print (C + (delta * K_baseline) - prodfunc(K_baseline, L_dict, global_dict))
print ('#############################################')

############## Counterfacutal of no immigrants ###############################

pop = {}

pop.update({'PH': 0.2})
pop.update({'RH': 0.005})
pop.update({'AH': 0.505})
pop.update({'PC': 0.02})
pop.update({'RC': 0.05})
pop.update({'AC': 0.32})
global_dict.update({'pop': pop})

HH_res_noimm, K_noimm, L_dict = OLG_solve_steady(global_dict)

w_noimm , r_noimm = prices(K_noimm, L_dict, global_dict)

print ('wages')
print ('wages test',w_noimm)
print ('wages baseline',w_baseline)
print ('interest rate')
print ('r test',r_noimm)
print ('r baseline',r_baseline)
print ('#############################################')
print ('noimm complete')
print ('HH_res',HH_res_noimm)
print ('a_max',K_noimm * gen_par['a grid scale'])
print ('#############################################')

# util_noimm = {}
# util_noimm['AH'] = np.sum(util(HH_res['AH'][2], HH_res['AH'][0], global_dict))
# util_noimm['AC'] = np.sum(util(HH_res['AC'][2], HH_res['AC'][0], global_dict))

L_noimm = dict.copy(L_dict)

Y_noimm = prodfunc(K,L_noimm, global_dict)
L_agg_noimm = agg_lab(L_noimm, global_dict)

print (HH_res_noimm)
print (w_noimm)
print (w_baseline)

quit()
print ('################## No immigrants ########################################')
print ('Wage percent change College',(w_noimm['AC'] - w_baseline['AC']) / w_baseline['AC'])
print ('Wage percent change HS',(w_noimm['AH'] - w_baseline['AH']) / w_baseline['AH'])
print ('Interest rate change', (r_noimm - r_baseline))
print ('Welfare change college', (util_noimm['AC'] - util_baseline['AC']) / np.abs(util_baseline['AC']))
print ('Welfare change HS', (util_noimm['AH'] - util_baseline['AH']) / np.abs(util_baseline['AH']))
print ('Labor change college', (L_noimm['AC'] - L_baseline['AC']) / L_baseline['AC'])
print ('Labor change HS', (L_noimm['AH'] - L_baseline['AH']) / L_baseline['AH'])
print ('Output change', (Y_noimm - Y_baseline) / Y_baseline)
print ('Agg labour change', (L_agg_noimm - L_agg_baseline) / L_agg_baseline)
print ('Agg capital change', (K - K_baseline) / K_baseline)
print ('############################################################################')

quit()

############### Counterfacutal of no immigrants ###############################

# pop = {}
#
# pop.update({'PH': 0.0})
# pop.update({'RH': 0.00})
# pop.update({'AH': 0.505})
# pop.update({'PC': 0.00})
# pop.update({'RC': 0.00})
# pop.update({'AC': 0.32})
#
# global_dict.update({'pop': pop})
#
# HH_res, K, L_dict = OLG_solve_steady(global_dict)
#
# w_noimm , r_noimm = prices(K, L_dict, global_dict)
#
# util_noimm = {}
# util_noimm['AH'] = np.sum(util(HH_res['AH'][2], HH_res['AH'][0], global_dict))
# util_noimm['AC'] = np.sum(util(HH_res['AC'][2], HH_res['AC'][0], global_dict))
#
# L_noimm = dict.copy(L_dict)
#
# Y_noimm = prodfunc(K,L_noimm, global_dict)
# L_agg_noimm = agg_lab(L_noimm, global_dict)
#
# print ('################## No immigrants ########################################')
# print ('Wage percent change College',(w_noimm['AC'] - w_baseline['AC']) / w_baseline['AC'])
# print ('Wage percent change HS',(w_noimm['AH'] - w_baseline['AH']) / w_baseline['AH'])
# print ('Interest rate change', (r_noimm - r_baseline))
# print ('Welfare change college', (util_noimm['AC'] - util_baseline['AC']) / np.abs(util_baseline['AC']))
# print ('Welfare change HS', (util_noimm['AH'] - util_baseline['AH']) / np.abs(util_baseline['AH']))
# print ('Labor change college', (L_noimm['AC'] - L_baseline['AC']) / L_baseline['AC'])
# print ('Labor change HS', (L_noimm['AH'] - L_baseline['AH']) / L_baseline['AH'])
# print ('Output change', (Y_noimm - Y_baseline) / Y_baseline)
# print ('Agg labour change', (L_agg_noimm - L_agg_baseline) / L_agg_baseline)
# print ('Agg capital change', (K - K_baseline) / K_baseline)
# print ('############################################################################')
#
# quit()

############### Counterfacutal of no immigrants same pop ###############################

# pop = {}
#
# pop.update({'PH': 0.00})
# pop.update({'RH': 0.00})
# pop.update({'AH': 0.62})
# pop.update({'PC': 0.00})
# pop.update({'RC': 0.00})
# pop.update({'AC': 0.38})
#
# global_dict.update({'pop': pop})
#
# HH_res, K, L_dict = OLG_solve_steady(global_dict)
#
# w_noimm , r_noimm = prices(K, L_dict, global_dict)
#
# util_noimm = {}
# util_noimm['AH'] = np.sum(util(HH_res['AH'][2], HH_res['AH'][0], global_dict))
# util_noimm['AC'] = np.sum(util(HH_res['AC'][2], HH_res['AC'][0], global_dict))
#
# L_noimm = dict.copy(L_dict)
#
# Y_noimm = prodfunc(K,L_noimm, global_dict)
# L_agg_noimm = agg_lab(L_noimm, global_dict)
#
# print ('################## No immigrants same pop ########################################')
# print ('Wage percent change College',(w_noimm['AC'] - w_baseline['AC']) / w_baseline['AC'])
# print ('Wage percent change HS',(w_noimm['AH'] - w_baseline['AH']) / w_baseline['AH'])
# print ('Interest rate change', (r_noimm - r_baseline))
# print ('Welfare change college', (util_noimm['AC'] - util_baseline['AC']) / np.abs(util_baseline['AC']))
# print ('Welfare change HS', (util_noimm['AH'] - util_baseline['AH']) / np.abs(util_baseline['AH']))
# print ('Labor change college', (L_noimm['AC'] - L_baseline['AC']) / L_baseline['AC'])
# print ('Labor change HS', (L_noimm['AH'] - L_baseline['AH']) / L_baseline['AH'])
# print ('Output change', (Y_noimm - Y_baseline) / Y_baseline)
# print ('Agg labour change', (L_agg_noimm - L_agg_baseline) / L_agg_baseline)
# print ('Agg capital change', (K - K_baseline) / K_baseline)
# print ('############################################################################')
#
# quit()

################ Counterfacutal of only Skilled immigrants ###############################
# pop = {}
#
# pop.update({'PH': 0.0})
# pop.update({'RH': 0.00})
# pop.update({'AH': 0.505})
# pop.update({'PC': 0.00})
# pop.update({'RC': 0.05})
# pop.update({'AC': 0.32})
#
# global_dict.update({'pop': pop})
#
# HH_res, K, L_dict = OLG_solve_steady(global_dict)
#
# w_Sonly , r_Sonly = prices(K, L_dict, global_dict)
#
# util_Sonly = {}
# util_Sonly['AH'] = np.sum(util(HH_res['AH'][2], HH_res['AH'][0], global_dict))
# util_Sonly['AC'] = np.sum(util(HH_res['AC'][2], HH_res['AC'][0], global_dict))
#
# L_Sonly = dict.copy(L_dict)
#
# Y_Sonly = prodfunc(K,L_Sonly, global_dict)
# L_agg_Sonly = agg_lab(L_Sonly, global_dict)
#
# print ('################## only Skilled immigrants ########################################')
# print ('Wage percent change College',(w_Sonly['AC'] - w_baseline['AC']) / w_baseline['AC'])
# print ('Wage percent change HS',(w_Sonly['AH'] - w_baseline['AH']) / w_baseline['AH'])
# print ('Interest rate change', (r_Sonly - r_baseline))
# print ('Welfare change college', (util_Sonly['AC'] - util_baseline['AC']) / np.abs(util_baseline['AC']))
# print ('Welfare change HS', (util_Sonly['AH'] - util_baseline['AH']) / np.abs(util_baseline['AH']))
# print ('Labor change college', (L_Sonly['AC'] - L_baseline['AC']) / L_baseline['AC'])
# print ('Labor change HS', (L_Sonly['AH'] - L_baseline['AH']) / L_baseline['AH'])
# print ('Output change', (Y_Sonly - Y_baseline) / Y_baseline)
# print ('Agg labour change', (L_agg_Sonly - L_agg_baseline) / L_agg_baseline)
# print ('Agg capital change', (K - K_baseline) / K_baseline)
# print ('############################################################################')
#
# quit()
# ############### Counterfacutal of anyone with a college degree ###############################
# pop = {}
#
# pop.update({'PH': 0.0})
# pop.update({'RH': 0.00})
# pop.update({'AH': 0.505})
# pop.update({'PC': 0.02})
# pop.update({'RC': 0.05})
# pop.update({'AC': 0.32})
#
# global_dict.update({'pop': pop})
#
# HH_res, K, L_dict = OLG_solve_steady(global_dict)
#
# w_coll , r_coll = prices(K, L_dict, global_dict)
#
# util_coll = {}
# util_coll['AH'] = np.sum(util(HH_res['AH'][2], HH_res['AH'][0], global_dict))
# util_coll['AC'] = np.sum(util(HH_res['AC'][2], HH_res['AC'][0], global_dict))
#
# L_coll = dict.copy(L_dict)
#
# Y_coll = prodfunc(K,L_coll, global_dict)
# L_agg_coll = agg_lab(L_coll, global_dict)
#
# print ('################## College ########################################')
# print ('Wage percent change College',(w_coll['AC'] - w_baseline['AC']) / w_baseline['AC'])
# print ('Wage percent change HS',(w_coll['AH'] - w_baseline['AH']) / w_baseline['AH'])
# print ('Interest rate change', (r_coll - r_baseline))
# print ('Welfare change college', (util_coll['AC'] - util_baseline['AC']) / np.abs(util_baseline['AC']))
# print ('Welfare change HS', (util_coll['AH'] - util_baseline['AH']) / np.abs(util_baseline['AH']))
# print ('Labor change college', (L_coll['AC'] - L_baseline['AC']) / L_baseline['AC'])
# print ('Labor change HS', (L_coll['AH'] - L_baseline['AH']) / L_baseline['AH'])
# print ('Output change', (Y_coll - Y_baseline) / Y_baseline)
# print ('Agg labour change', (L_agg_coll - L_agg_baseline) / L_agg_baseline)
# print ('Agg capital change', (K - K_baseline) / K_baseline)
# print ('############################################################################')

# ############### Only unskilled ###############################
pop = {}

pop.update({'PH': 0.105})
pop.update({'RH': 0.005})
pop.update({'AH': 0.505})
pop.update({'PC': 0.00})
pop.update({'RC': 0.00})
pop.update({'AC': 0.32})

global_dict.update({'pop': pop})

HH_res, K, L_dict = OLG_solve_steady(global_dict)

w_Uonly , r_Uonly = prices(K, L_dict, global_dict)

util_Uonly = {}
util_Uonly['AH'] = np.sum(util(HH_res['AH'][2], HH_res['AH'][0], global_dict))
util_Uonly['AC'] = np.sum(util(HH_res['AC'][2], HH_res['AC'][0], global_dict))

L_Uonly = dict.copy(L_dict)

Y_Uonly = prodfunc(K,L_Uonly, global_dict)
L_agg_Uonly = agg_lab(L_Uonly, global_dict)

print ('################## Unskilled only ########################################')
print ('Wage percent change Uonlyege',(w_Uonly['AC'] - w_baseline['AC']) / w_baseline['AC'])
print ('Wage percent change HS',(w_Uonly['AH'] - w_baseline['AH']) / w_baseline['AH'])
print ('Interest rate change', (r_Uonly - r_baseline))
print ('Welfare change Uonlyege', (util_Uonly['AC'] - util_baseline['AC']) / np.abs(util_baseline['AC']))
print ('Welfare change HS', (util_Uonly['AH'] - util_baseline['AH']) / np.abs(util_baseline['AH']))
print ('Labor change Uonlyege', (L_Uonly['AC'] - L_baseline['AC']) / L_baseline['AC'])
print ('Labor change HS', (L_Uonly['AH'] - L_baseline['AH']) / L_baseline['AH'])
print ('Output change', (Y_Uonly - Y_baseline) / Y_baseline)
print ('Agg labour change', (L_agg_Uonly - L_agg_baseline) / L_agg_baseline)
print ('Agg capital change', (K - K_baseline) / K_baseline)
print ('############################################################################')

quit()

################ Counterfacutal of only Skilled immigrants with same visas ###############################
pop = {}

pop.update({'PH': 0.0})
pop.update({'RH': 0.00})
pop.update({'AH': 0.505})
pop.update({'PC': 0.00})
pop.update({'RC': 0.155})
pop.update({'AC': 0.32})

global_dict.update({'pop': pop})

HH_res, K, L_dict = OLG_solve_steady(global_dict)

w_Sonly , r_Sonly = prices(K, L_dict, global_dict)

util_Sonly = {}
util_Sonly['AH'] = np.sum(util(HH_res['AH'][2], HH_res['AH'][0], global_dict))
util_Sonly['AC'] = np.sum(util(HH_res['AC'][2], HH_res['AC'][0], global_dict))

L_Sonly = dict.copy(L_dict)

Y_Sonly = prodfunc(K,L_Sonly, global_dict)
L_agg_Sonly = agg_lab(L_Sonly, global_dict)

print ('################## only Skilled immigrants same visas ########################################')
print ('Wage percent change College',(w_Sonly['AC'] - w_baseline['AC']) / w_baseline['AC'])
print ('Wage percent change HS',(w_Sonly['AH'] - w_baseline['AH']) / w_baseline['AH'])
print ('Interest rate change', (r_Sonly - r_baseline))
print ('Welfare change college', (util_Sonly['AC'] - util_baseline['AC']) / np.abs(util_baseline['AC']))
print ('Welfare change HS', (util_Sonly['AH'] - util_baseline['AH']) / np.abs(util_baseline['AH']))
print ('Labor change college', (L_Sonly['AC'] - L_baseline['AC']) / L_baseline['AC'])
print ('Labor change HS', (L_Sonly['AH'] - L_baseline['AH']) / L_baseline['AH'])
print ('Output change', (Y_Sonly - Y_baseline) / Y_baseline)
print ('Agg labour change', (L_agg_Sonly - L_agg_baseline) / L_agg_baseline)
print ('Agg capital change', (K - K_baseline) / K_baseline)
print ('############################################################################')

# ############### Counterfacutal of anyone with a college degree ###############################
pop = {}

pop.update({'PH': 0.0})
pop.update({'RH': 0.00})
pop.update({'AH': 0.505})
pop.update({'PC': 0.03})
pop.update({'RC': 0.125})
pop.update({'AC': 0.32})

global_dict.update({'pop': pop})

HH_res, K, L_dict = OLG_solve_steady(global_dict)

w_coll , r_coll = prices(K, L_dict, global_dict)

util_coll = {}
util_coll['AH'] = np.sum(util(HH_res['AH'][2], HH_res['AH'][0], global_dict))
util_coll['AC'] = np.sum(util(HH_res['AC'][2], HH_res['AC'][0], global_dict))

L_coll = dict.copy(L_dict)

Y_coll = prodfunc(K,L_coll, global_dict)
L_agg_coll = agg_lab(L_coll, global_dict)

print ('################## College same visas ########################################')
print ('Wage percent change College',(w_coll['AC'] - w_baseline['AC']) / w_baseline['AC'])
print ('Wage percent change HS',(w_coll['AH'] - w_baseline['AH']) / w_baseline['AH'])
print ('Interest rate change', (r_coll - r_baseline))
print ('Welfare change college', (util_coll['AC'] - util_baseline['AC']) / np.abs(util_baseline['AC']))
print ('Welfare change HS', (util_coll['AH'] - util_baseline['AH']) / np.abs(util_baseline['AH']))
print ('Labor change college', (L_coll['AC'] - L_baseline['AC']) / L_baseline['AC'])
print ('Labor change HS', (L_coll['AH'] - L_baseline['AH']) / L_baseline['AH'])
print ('Output change', (Y_coll - Y_baseline) / Y_baseline)
print ('Agg labour change', (L_agg_coll - L_agg_baseline) / L_agg_baseline)
print ('Agg capital change', (K - K_baseline) / K_baseline)
print ('############################################################################')

# ############### Only unskilled same visas ###############################
pop = {}

pop.update({'PH': 0.17})
pop.update({'RH': 0.01})
pop.update({'AH': 0.505})
pop.update({'PC': 0.00})
pop.update({'RC': 0.00})
pop.update({'AC': 0.32})

global_dict.update({'pop': pop})

HH_res, K, L_dict = OLG_solve_steady(global_dict)

w_Uonly , r_Uonly = prices(K, L_dict, global_dict)

util_Uonly = {}
util_Uonly['AH'] = np.sum(util(HH_res['AH'][2], HH_res['AH'][0], global_dict))
util_Uonly['AC'] = np.sum(util(HH_res['AC'][2], HH_res['AC'][0], global_dict))

L_Uonly = dict.copy(L_dict)

Y_Uonly = prodfunc(K,L_Uonly, global_dict)
L_agg_Uonly = agg_lab(L_Uonly, global_dict)

print ('################## Unskilled only same visas ########################################')
print ('Wage percent change Uonlyege',(w_Uonly['AC'] - w_baseline['AC']) / w_baseline['AC'])
print ('Wage percent change HS',(w_Uonly['AH'] - w_baseline['AH']) / w_baseline['AH'])
print ('Interest rate change', (r_Uonly - r_baseline))
print ('Welfare change Uonlyege', (util_Uonly['AC'] - util_baseline['AC']) / np.abs(util_baseline['AC']))
print ('Welfare change HS', (util_Uonly['AH'] - util_baseline['AH']) / np.abs(util_baseline['AH']))
print ('Labor change Uonlyege', (L_Uonly['AC'] - L_baseline['AC']) / L_baseline['AC'])
print ('Labor change HS', (L_Uonly['AH'] - L_baseline['AH']) / L_baseline['AH'])
print ('Output change', (Y_Uonly - Y_baseline) / Y_baseline)
print ('Agg labour change', (L_agg_Uonly - L_agg_baseline) / L_agg_baseline)
print ('Agg capital change', (K - K_baseline) / K_baseline)
print ('############################################################################')
