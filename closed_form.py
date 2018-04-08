import numpy as np
import scipy
from scipy import optimize

life = 2
na = 11

# Utility function is the same storelsatten

util_par = {}
util_par.update({'beta': 1.0})
util_par.update({'alpha': 1.0})
util_par.update({'gamma': 4.0})

# Psi is the proportion of the population for each group
psi = {}
psi.update({'S': 1.0})

prod_par = {}
prod_par.update({'A': 1.0})
prod_par.update({'theta': 0.4})
prod_par.update({'delta': (((prod_par['theta'] * (1.0 / (3.0 / 5.0))) - ((1.05 ** 5) - 1)))})

pop_generic = (1.0 / life) * np.ones(life)

pop = {}
pop.update({'S': psi['S'] * pop_generic})

eff = {}
eff.update({'S': 1.0})

types = len(list(pop.keys()))
types_key = list(pop.keys())

w = 1.0
r = 0.5

global_dict = {}
global_dict.update({'prod_par': prod_par})
global_dict.update({'util_par': util_par})
global_dict.update({'psi': psi})
global_dict.update({'pop': pop})
global_dict.update({'eff': eff})
global_dict.update({'life': life})
global_dict.update({'types_key': types_key})
global_dict.update({'points': na})
global_dict.update({'r':r})
global_dict.update({'w':w})

##### closed form solution for two periods ###########

def c(a_0, global_dict):

    util_par = global_dict['util_par']
    eff = global_dict['eff']
    w = global_dict['w']
    r = global_dict['r']

    tau_l = global_dict['l_tax']
    tau_k = global_dict['k_tax']

    cons1 = 1 + util_par['beta'] + util_par['alpha'] + (util_par['beta'] * util_par['alpha'])

    c_0 = (1.0 / cons1) * (((w * eff['S'] * (1 - tau_l)) * (1.0 + (1.0 / (1.0 + (r * (1 - tau_k)))))) + ((1.0 + r) * a_0))
    c_1 = util_par['beta'] * (1 + (r * (1 - tau_k))) * c_0

    return c_0, c_1

def l(a_0, global_dict):

    util_par = global_dict['util_par']
    eff = global_dict['eff']
    w = global_dict['w']
    r = global_dict['r']

    tau_l = global_dict['l_tax']
    tau_k = global_dict['k_tax']

    l_0 = 1 - ((util_par['alpha'] * c(a_0, global_dict)[0]) / (w * eff['S'] * (1 - tau_l)))
    l_1 = 1 - ((util_par['alpha'] * c(a_0, global_dict)[1]) / (w * eff['S'] * (1 - tau_l)))

    return l_0, l_1

def a_p(a_0, global_dict):

    eff = global_dict['eff']
    w = global_dict['w']
    r = global_dict['r']

    tau_l = global_dict['l_tax']
    tau_k = global_dict['k_tax']

    a_1 = (1.0 / (1.0 + (r * (1 - tau_k)))) * (c(a_0, global_dict)[1] - (w * (1 - tau_l) * eff['S'] * l(a_0, global_dict)[1]))

    return a_1

def solve(a_0, global_dict):

    a_1 = a_p(a_0, global_dict)

    c_seq = c(a_0, global_dict)
    l_seq = l(a_0, global_dict)

    return a_0, l_seq, a_1, c_seq

global_dict.update({'l_tax' : 0.2})
global_dict.update({'k_tax' : 0.2})

print (solve(0.0, global_dict))

##### closed form solution for two periods of working with bequests ###########
#
# def c(a_0, global_dict):
#
#     util_par = global_dict['util_par']
#     eff = global_dict['eff']
#     w = global_dict['w']
#     r = global_dict['r']
#     beq = global_dict['beq']
#
#     cons1 = 1 + util_par['beta'] + util_par['alpha'] + (util_par['beta'] * util_par['alpha'])
#
#     c_0 = (1.0 / cons1) * (((w * eff['S']) * (1.0 + (1.0 / ((1.0 + r) * (1.0 - beq))))) + ((1.0 + r) * a_0))
#     c_1 = util_par['beta'] * (1.0 + r) * (1.0 - beq) * c_0
#
#     return c_0, c_1
#
# def l(a_0, global_dict):
#
#     util_par = global_dict['util_par']
#     eff = global_dict['eff']
#     w = global_dict['w']
#     r = global_dict['r']
#
#     l_0 = 1 - ((util_par['alpha'] * c(a_0, global_dict)[0]) / (w * eff['S']))
#     l_1 = 1 - ((util_par['alpha'] * c(a_0, global_dict)[1]) / (w * eff['S']))
#
#     return l_0, l_1
#
# def a_p(a_0, global_dict):
#
#     eff = global_dict['eff']
#     w = global_dict['w']
#     r = global_dict['r']
#     beq = global_dict['beq']
#
#     a_1conso = (1.0 / ((1.0 + r) * (1 - beq)))
#     a_1consi = c(a_0, global_dict)[1] - ((w * eff['S']) * l(a_0, global_dict)[1])
#
#     a_1 = a_1conso * a_1consi
#
#     return a_1
#
# def solve(a_0, global_dict):
#
#     a_1 = a_p(a_0, global_dict)
#
#     c_seq = c(a_0, global_dict)
#     l_seq = l(a_0, global_dict)
#
#     return a_0, l_seq, c_seq, a_1
#
# beq = 0.2
# global_dict.update({'beq': beq})
#
# a1s_conso1 = (util_par['beta'] * (1.0 + util_par['alpha']) * beq * ((1.0 + r) ** 2))
# a1s_conso2 = (1.0 + util_par['beta']) * (1 + util_par['alpha'])
#
# a1s_conso = 1 - (a1s_conso1 / a1s_conso2)
#
# a1s_consi1 = (util_par['beta'] * (1.0 + util_par['alpha'])) / ((1.0 + util_par['beta']) * (1 + util_par['alpha']))
# a1s_consi2 = 1.0 / ((1.0 + r) * (1.0 - beq))
#
# a1s_consi = (a1s_consi1 * (1 + a1s_consi2)) - a1s_consi2
#
# a1s = (1.0 / a1s_conso) * (w * eff['S']) * a1s_consi
#
# print (a1s)
#
# print (solve(beq * (1.0 + r) * a1s, global_dict))

##### closed form solution for two periods ###########
#
# def c(a_0, global_dict):
#
#     util_par = global_dict['util_par']
#     eff = global_dict['eff']
#     w = global_dict['w']
#     r = global_dict['r']
#
#     cons1 = 1 + util_par['beta'] + util_par['alpha'] + (util_par['beta'] * util_par['alpha'])
#
#     c_0 = (1.0 / cons1) * (((w * eff['S']) * (1.0 + (1.0 / (1.0 + r)))) + ((1.0 + r) * a_0))
#     c_1 = util_par['beta'] * (1 + r) * c_0
#
#     return c_0, c_1
#
# def l(a_0, global_dict):
#
#     util_par = global_dict['util_par']
#     eff = global_dict['eff']
#     w = global_dict['w']
#     r = global_dict['r']
#
#     l_0 = 1 - ((util_par['alpha'] * c(a_0, global_dict)[0]) / (w * eff['S']))
#     l_1 = 1 - ((util_par['alpha'] * c(a_0, global_dict)[1]) / (w * eff['S']))
#
#     return l_0, l_1
#
# def a_p(a_0, global_dict):
#
#     eff = global_dict['eff']
#     w = global_dict['w']
#     r = global_dict['r']
#
#     a_1 = (1.0 / (1.0 + r)) * (c(a_0, global_dict)[1] - (w * eff['S'] * l(a_0, global_dict)[1]))
#
#     return a_1
#
# def solve(a_0, global_dict):
#
#     a_1 = a_p(a_0, global_dict)
#
#     c_seq = c(a_0, global_dict)
#     l_seq = l(a_0, global_dict)
#
#     return a_0, l_seq, a_1, c_seq
#
# print (c(0.0, global_dict))
# cons1 = 1 + util_par['beta'] + util_par['alpha'] + (util_par['beta'] * util_par['alpha'])
# print (cons1)
# print c_0 = (1.0 / cons1) * (((w * eff['S']) * (1.0 + (1.0 / (1.0 + r)))) + ((1.0 + r) * a_0))

###### closed form solution for three periods tw in work one in retirement ###########

# def c(a_0, global_dict):
#
#     util_par = global_dict['util_par']
#     eff = global_dict['eff']
#     w = global_dict['w']
#     r = global_dict['r']
#
#     cons1 = 1 + util_par['beta'] + (util_par['beta'] ** 2) + util_par['alpha'] + (util_par['alpha'] * util_par['beta'])
#
#     c_0 = (1.0 / cons1) * (((w * eff['S']) * (1.0 + (1.0 / (1.0 + r)))) + ((1.0 + r) * a_0))
#     c_1 = util_par['beta'] * (1 + r) * c_0
#     c_2 = util_par['beta'] * (1 + r) * c_1
#
#     return c_0, c_1, c_2
#
# def l(a_0, global_dict):
#
#     util_par = global_dict['util_par']
#     eff = global_dict['eff']
#     w = global_dict['w']
#     r = global_dict['r']
#
#     l_0 = 1 - ((util_par['alpha'] * c(a_0, global_dict)[0]) / (w * eff['S']))
#     l_1 = 1 - ((util_par['alpha'] * c(a_0, global_dict)[1]) / (w * eff['S']))
#     l_2 = 0
#
#     return l_0, l_1, l_2
#
# def a_p(a_0, global_dict):
#
#     eff = global_dict['eff']
#     w = global_dict['w']
#     r = global_dict['r']
#
#     a_1 = (w * eff['S'] * l(a_0, global_dict)[0]) + ((1 + r) * a_0) - c(a_0, global_dict)[0]
#     a_2 = c(a_0, global_dict)[2] / (1.0 + r)
#
#     return a_1, a_2
#
# def solve(a_0, global_dict):
#
#     a_seq = a_p(a_0, global_dict)
#
#     c_seq = c(a_0, global_dict)
#     l_seq = l(a_0, global_dict)
#
#     return a_0, l_seq, a_seq, c_seq
#
# print (solve(0.0, global_dict))
#

# ###### closed form solution for two periods one in retirement with bequests ###########

# def c(a_0, global_dict):
#
#     util_par = global_dict['util_par']
#     eff = global_dict['eff']
#     w = global_dict['w']
#     r = global_dict['r']
#     beq = global_dict['b']
#
#     cons1 = 1 + util_par['beta'] + util_par['alpha']
#
#     c_0 = (1.0 / cons1) * ((w * eff['S']) + ((1.0 + r) * a_0))
#     c_1 = util_par['beta'] * (1 + r) * (1 - beq) * c_0
#
#     return c_0, c_1
#
# def l(a_0, global_dict):
#
#     util_par = global_dict['util_par']
#     eff = global_dict['eff']
#     w = global_dict['w']
#     r = global_dict['r']
#
#     l_0 = 1 - ((util_par['alpha'] * c(a_0, global_dict)[0]) / (w * eff['S']))
#     l_1 = 0
#
#     return l_0, l_1
#
# def a_p(a_0, global_dict):
#
#     eff = global_dict['eff']
#     w = global_dict['w']
#     r = global_dict['r']
#     beq = global_dict['b']
#
#     a_1 = (1.0 / ((1 + r) * (1 - beq))) * c(a_0, global_dict)[0]
#
#     return a_1
#
# def solve(a_0, global_dict):
#
#     a_1 = a_p(a_0, global_dict)
#
#     c_seq = c(a_0, global_dict)
#     l_seq = l(a_0, global_dict)
#
#     return a_0, l_seq, c_seq, a_1
#
# # print (solve(0.0, global_dict))
#
# beq = 0.2
#
# global_dict.update({'b':beq})
#
# a1s_den = util_par['beta'] * w * eff['S']
# a1s_num = (1 + util_par['beta'] + util_par['alpha']) - (util_par['beta'] * beq * ((1 + r) ** 2))
#
# a_1_star = a1s_den / a1s_num
#
# print (solve((1 + r) * beq * a_1_star, global_dict))
#
# quit()
###### closed form solution for three periods ###########

# def c(a_0, global_dict):
#
#     util_par = global_dict['util_par']
#     eff = global_dict['eff']
#     w = global_dict['w']
#     r = global_dict['r']
#
#     cons1 = 1 + util_par['beta'] + (util_par['beta'] ** 2) \
#                 + util_par['alpha'] + (util_par['beta'] * util_par['alpha']) + ((util_par['beta'] ** 2) * util_par['alpha'])
#
#     consr = (1.0 / (1.0 + r))
#     consw = (w * eff['S']) * (1 + consr + (consr ** 2))
#
#     c_0 = (1.0 / cons1) * (consw + ((1.0 + r) * a_0))
#     c_1 = util_par['beta'] * (1 + r) * c_0
#     c_2 = util_par['beta'] * (1 + r) * c_1
#
#     return c_0, c_1, c_2
#
# def l(a_0, global_dict):
#
#     util_par = global_dict['util_par']
#     eff = global_dict['eff']
#     w = global_dict['w']
#     r = global_dict['r']
#
#     l_0 = 1 - ((util_par['alpha'] * c(a_0, global_dict)[0]) / (w * eff['S']))
#     l_1 = 1 - ((util_par['alpha'] * c(a_0, global_dict)[1]) / (w * eff['S']))
#     l_2 = 1 - ((util_par['alpha'] * c(a_0, global_dict)[2]) / (w * eff['S']))
#
#     return l_0, l_1, l_2
#
# def a_p(a_0, global_dict):
#
#     eff = global_dict['eff']
#     w = global_dict['w']
#     r = global_dict['r']
#
#     a_1 = (1.0 / (1.0 + r)) * (c(a_0, global_dict)[1] - (w * eff['S'] * l(a_0, global_dict)[1])) \
#         + ((1.0 / (1.0 + r)) ** 2) * (c(a_0, global_dict)[2] - (w * eff['S'] * l(a_0, global_dict)[2]))
#
#     a_2 = (1.0 / (1.0 + r)) * (c(a_0, global_dict)[2] - (w * eff['S'] * l(a_0, global_dict)[2]))
#
#     return a_1, a_2
#
# def solve(a_0, global_dict):
#
#     a_seq = a_p(a_0, global_dict)
#
#     c_seq = c(a_0, global_dict)
#     l_seq = l(a_0, global_dict)
#
#     return a_0, l_seq, a_seq, c_seq
#
# print (solve(0.0, global_dict))

##### as a test I find the intial assets that mean you dont want to work in firrst period ########
# icons1 = (1 + util_par['alpha']) * (1 + util_par['beta'] + (util_par['beta'] ** 2))
# consr = (1.0 / (1.0 + r))
# icons2 = 1 + consr + (consr ** 2)
#
# a_0_star = (w * eff['S'] / (1.0 + r)) * (((1.0 / util_par['alpha']) * icons1) - icons2)
#
# print (a_0_star)
#
# print (l(a_0_star, global_dict))
