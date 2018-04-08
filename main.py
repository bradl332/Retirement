import numpy as np
import scipy
from scipy import optimize
from scipy.interpolate import make_interp_spline
import multiprocessing as mp
import pickle

from OLG_solve_steady import OLG_solve_steady
from calibration import calibration_util, calibration_eps
from funcs import prodfunc, prices

work = 9
retire = 5
life = work + retire
na = 500

util_par = {}
util_par.update({'alpha': 1.779})
util_par.update({'beta': 0.835})

rho = {}
rho.update({'SU': 10000000.0})

lamda = {}
lamda_su = {}
lamda_u = {}
lamda_s = {}

lamda_su.update({'S': 1.0})
lamda_su.update({'U': 1.0})
lamda_su.update({'M': 1.0})

lamda.update({'SU': lamda_su})

prod_par = {}
prod_par.update({'theta': 0.36})

per_length = 70.0 / life
delta = (prod_par['theta'] * (1.0 / (3.0 / per_length)) - (((1.05) ** (per_length)) - 1))

prod_par.update({'delta': delta})
prod_par.update({'A': 1.0})
prod_par.update({'lamda': lamda})
prod_par.update({'rho': rho})

gen_par = {}
gen_par.update({'bequest': 0.2})
gen_par.update({'a grid scale': 4.0})
gen_par.update({'pop_growth': 0.01})

eff = {}

eff.update({'S': 2.0})
eff.update({'U': 1.0})
eff.update({'M': 1.5})

pop = {}

pop.update({'S': 0.33})
pop.update({'U': 0.33})
pop.update({'M': 0.33})

pop_age = {}
pop_age.update({'S': pop['S'] * (1.0 / life) * np.ones(life)})
pop_age.update({'U': pop['U'] * (1.0 / life) * np.ones(life)})
pop_age.update({'M': pop['M'] * (1.0 / life) * np.ones(life)})

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

HH_res, K, L_dict = OLG_solve_steady(global_dict)

print (K / prodfunc(K, L_dict, global_dict))
print ('HH_res baseline',HH_res)
print ('a_max baseline',global_dict['gen_par']['a grid scale'] * K)
print (prices(K, L_dict, global_dict))
print('#####################')

quit()
# w,r = prices(K, L_dict, global_dict)
# print ('K_Y',K / prodfunc(K, L_dict, global_dict))
# print ('Average hours PH',np.dot(pop_age['PH'], HH_res['PH'][0]))
# print ('Average hours RH',np.dot(pop_age['RH'], HH_res['RH'][0]))
# print ('Average hours AH',np.dot(pop_age['AH'], HH_res['AH'][0]))
# print ('Average hours PC',np.dot(pop_age['PC'], HH_res['PC'][0]))
# print ('Average hours RC',np.dot(pop_age['RC'], HH_res['RC'][0]))
# print ('Average hours AC',np.dot(pop_age['AC'], HH_res['AC'][0]))
# print ('rental rate',r)
# print ("")
# # print ((w['AH'] * eff['AH']) * np.dot(pop_age['AH'], HH_res['AH'][0]))
# # print (0.8 * prodfunc(K, L_dict, global_dict))
# print ("")
# print ('AH earnings target 1.37', (w['AH'] * eff['AH']) / (w['PH'] * eff['PH']))
# print ('RH earnings target 1.48', (w['RH'] * eff['RH']) / (w['PH'] * eff['PH']))
# print ('PC earnings target 1.48', (w['PC'] * eff['PC']) / (w['PH'] * eff['PH']))
# print ('RC earnings target 2.67', (w['RC'] * eff['RC']) / (w['PH'] * eff['PH']))
# print ('AC earnings target 2,24', (w['AC'] * eff['AC']) / (w['PH'] * eff['PH']))
#
# quit()

# I calibrate beta and alpha on their own #

# calib_util = scipy.optimize.fsolve(calibration_util, np.array([0.925, 1.88]), global_dict)

# Now to calibrate the epsilons

# bounds = []
#
# bounds.append((0.85,0.995))
# bounds.append((1.5,2.2))
# bounds.append((1.0, 10.0))
# bounds.append((1.0, 10.0))
# bounds.append((1.0, 10.0))
# bounds.append((1.0, 10.0))
# bounds.append((1.0, 10.0))

# calib_eps = scipy.optimize.differential_evolution(calibration_eps, bounds, [global_dict])

calib_eps = scipy.optimize.fsolve(calibration_eps, np.asarray([ 0.85,  1.75,  1.97,  1.62,  1.68,  3.5, 2.9, 2.5, 2.14]), global_dict)
