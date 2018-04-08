import numpy as np
import scipy
import sys
from scipy.interpolate import make_interp_spline

from funcs import prodfunc, prices, util, agg_lab

def agent_solve(a_p, a, args):

    # This function takes current period a as given and calcualtes the optimal
    # choice of next period assets.

    tipo = args['types_key']

    w = args['wages']
    r = args['r']

    util_par = args['util_par']
    eff = args['eff']

    bq = args['bq']

    tau_l = args['lab tax']
    tau_k = args['cap tax']
    tra = args['tra']

    # This solves the retirement problem. It needs to be separate so I can force
    # the labour supply to 0.

    # Define Non labour income it makes adding in transfers etc much easier later on

    nlinc = ((1.0 + (r * (1 - tau_k))) * (a * bq)) - a_p

    if args['version'] == 'retire':

        c = nlinc
        l = 0.0 * c

    # using the intratemporal constraint I can solve for l. Again since a_p is
    # a vector we see all possible values of l for each choice of a_p. Some of these
    # values will violate labour being >= or <= 1. The code can make exceptions in
    # this case.

    else:

        #################################################################################################################################################
        # WARNING IF YOU CHANGE UTILITY FUNCTION (log c + alpha * log (1 - l)) YOU MUST CHANGE THIS
        # LINE OF CODE AND THE ASSOCIATED UTILITY FUNCTION

        conso = w[tipo] * eff[tipo] * (1.0 - tau_l) * (1.0 + util_par['alpha'])
        conso = (1.0 / conso)

        consi = (w[tipo] * eff[tipo] * (1.0 - tau_l)) - (util_par['alpha'] * (nlinc))

        l = conso * consi
        #################################################################################################################################################

        c = np.zeros(l.shape)

        # Here is where the code handles exceptions

        c[np.where(l <= 0)] = nlinc[np.where(l <= 0)]
        l[np.where(l <= 0)] = 0

        c[np.where(l >= 1)] = (w[tipo] * eff[tipo] * (1.0 - tau_l)) + (nlinc[np.where(l >= 1)])
        l[np.where(l >= 1)] = 1.0

        #################################################################################################################################################
        # WARNING IF YOU CHANGE UTILITY FUNCTION (log c + alpha * log (1 - l)) YOU MUST CHANGE THIS
        # LINE OF CODE AND THE ASSOCIATED UTILITY FUNCTION

        c[np.logical_and(l > 0.0, l < 1.0)] = (w[tipo] * eff[tipo] * (1.0 - tau_l)) * (1.0 / util_par['alpha']) * (1 - l[np.logical_and(l > 0.0, l < 1.0)])

        #################################################################################################################################################

    utility = util(c, l, args)

    # This handles the situation in case of 0 consumption (sometimes happens in the retirement calcualtion if the agent has 0 assets)

    utility[np.where(c < 0.0)] = - np.inf
    c[np.where(c < 0.0)] = 0.0

    # Create a vector of todays utility + continuation value. Remember c should be strictly decreasing
    # and V_nex should be strictly increasing.

    dec_vec = utility + (util_par['beta'] * args['V_nex'])

    # Tells us the location in the vector of the optimal choice of a_p

    a_p_loc = np.argmax(dec_vec)

    a_p = args['a_grid'][a_p_loc]
    V_n = dec_vec[a_p_loc]
    c = c[a_p_loc]
    l = l[a_p_loc]

    return a_p, V_n, c, l

def agent_simul(a_p, a, args):

    # This is very similar to agent_solve and could probably be nested within it.
    # I'm not a fan with how python handles stuff with technically 0 dimensions.

    tipo = args['types_key']

    w = args['wages']
    r = args['r']

    util_par = args['util_par']
    eff = args['eff']

    bq = args['bq']

    tau_l = args['lab tax']
    tau_k = args['cap tax']
    tra = args['tra']

    nlinc = ((1.0 + (r * (1 - tau_k))) * (a * bq)) - a_p

    if args['version'] == 'retire':

        c = nlinc
        l = 0.0 * c

    else:

        #################################################################################################################################################
        # WARNING IF YOU CHANGE UTILITY FUNCTION (log c + alpha * log (1 - l)) YOU MUST CHANGE THIS
        # LINE OF CODE AND THE ASSOCIATED UTILITY FUNCTION

        conso = w[tipo] * eff[tipo] * (1.0 - tau_l) * (1.0 + util_par['alpha'])
        conso = (1.0 / conso)

        consi = (w[tipo] * eff[tipo] * (1.0 - tau_l)) - (util_par['alpha'] * nlinc)

        l = conso * consi

        #################################################################################################################################################
        if l <= 0.0:
            l = 0.0
            c = nlinc

        if l >= 1.0:
            l = 1.0
            c = (w[tipo] * eff[tipo] * (1.0 - tau_l)) + (nlinc)

        #################################################################################################################################################
        # WARNING IF YOU CHANGE UTILITY FUNCTION (log c + alpha * log (1 - l)) YOU MUST CHANGE THIS
        # LINE OF CODE AND THE ASSOCIATED UTILITY FUNCTION

        if l < 1.0 and l > 0.0:
            c = ((w[tipo] * eff[tipo] * (1.0 - tau_l)) * (1.0 / util_par['alpha'])) * (1 - l)
        #################################################################################################################################################

    return c,l

def para_agent(y):

    # this function could be done away with. It stays in case I want to parallelize
    # more stuff. i refers to which point in the asset grid I want to solve the problem for.
    # This function takes the entire asset grid to calculate the value over
    # all future choices of a_p, a current period asset holding and all the
    # relevant arguments. It returns teh optimal choice of a_p, its associated
    # value. it also returns the optimal choice of c, l. I don't actually use
    # these they are there more in case I want to check stuff later on.

    i = y[0]
    arg = y[1]

    a_p, V_n, c, l = agent_solve(arg['a_grid'], arg['a_grid'][i], arg)

    return a_p, V_n, c, l

def para_type_steady(y):

    # Without the folloiwng line python gets pissy

    # y = y[0]

    na = y['points']
    life = y['life']
    work = y['work']
    retire = y['retire']
    tipo = y['types_key']

    if y['pop'][tipo] == 0.0:
        return np.zeros(life), np.zeros(life), np.zeros(life), tipo

    else:

        a_pol_list = []
        c_pol_list = []
        l_pol_list = []

        a_grid = y['a_grid']

        # Technically I could take a bequest each period, but that'd be daft
        # Instead I just limit to the end of life, to do this I say 0 bequests until
        # the last period ie you keep all your wealth, in the final period you lose
        # bequest amount of your wealth to be given to the new generation of your type

        bq = np.ones(life)
        bq[-1] = bq[-1] - y['gen_par']['bequest']
        y.update({'bq': - 5.0})

        y.update({'lab tax': None})
        y.update({'cap tax': None})
        y.update({'tra': None})

        for j in range(life - 1, - 1, -1):

            # The agent problem is set up differnetly for retirement and working
            # I have to force labour = 0 and that agents can only consume wealth
            # in retirement. In work they can do whatever they want.

            if j >= work:
                y['version'] = 'retire'
            else:
                y['version'] = 'work'

            para_args = []
            results = []

            # set the amount that is to be bequested in each period
            y['bq'] = bq[j]

            # Set the taxes and transfers for each period
            y['lab tax'] = y['taxes']['labor tax'][j]
            y['cap tax'] = y['taxes']['capital tax'][j]
            y['tra'] = y['taxes']['transfers'][j]

            # for each point possible a that the agent's can enter the state with
            # solve the agent's problem. This also works for the last period guys
            # since they have a V_nex and therefore will always choose a_p = 0
            # This again is set up in case I ever wanted to parallel this section
            # rather than over the types.

            for i in range(na):
                para_args.append((i,y))
                results.append(para_agent(para_args[i]))

            # results is returned as a list of tuples. The list(zip) command simply
            # converts it 4 containing, a_p, continuation value, consumption and labour choices
            results = list(zip(*results))

            a_pol = np.asarray(results[0])
            V_p = np.asarray(results[1])
            c_pol = np.asarray(results[2])
            l_pol = np.asarray(results[3])

            # Add the policy functions to the list. Each list will contain the policy
            # function at each age j. I could store V_nex aswell but doesn't seem much points

            a_pol_list.append(a_pol)
            c_pol_list.append(c_pol)
            l_pol_list.append(l_pol)

            V_nex = np.copy(V_p)

            y['V_nex'] = V_nex

        # Remember since this code solves the agent's problem backwards but the list
        # containing the policy functions has them in sequential order. The next lines reverse
        # the list.

        a_pol_list = a_pol_list[::-1]
        c_pol_list = c_pol_list[::-1]
        l_pol_list = l_pol_list[::-1]

        # Now I have the policy functions at each age I can give the agent an initial
        # asset holding and see what they will do. Given that c,l and are functions
        # of asset holdings it seems sensible to create the path of asset holdings
        # then using intratemporal condition and BC calcualte c,l. I use splines so that
        # I can start the a_path from any point within the range of the a_grid without
        # having to fart around having to adjust stuff.

        a_path = np.zeros(life)

        a_path[0] = 0.0

        for j in range(life - 1):
            a_pol_func = make_interp_spline(a_grid, a_pol_list[j])
            a_path[j + 1] = a_pol_func(a_path[j])

        # Given that bequests can be positive and passed to the youngest generation of
        # each type that will affect their choices of a_p which will affect how much
        # will be bequested. These next lines find the initial asset holdings that
        # will be consistent with the amount agents will be holding in the last period.

        bequest = (y['pop_age'][tipo][-1] / ((1.0 + y['gen_par']['pop_growth']) * y['pop_age'][tipo][0])) * (1 - bq[-1]) * (1 + y['r']) * a_path[-1]

        death_dist = 1.0
        death_tol = 0.0001
        death_lam = 0.8

        while death_dist > death_tol:
            a_path = np.zeros(life)
            a_path[0] = bequest

            for j in range(life - 1):
                a_pol_func = make_interp_spline(a_grid, a_pol_list[j])
                a_path[j + 1] = a_pol_func(a_path[j])

            bequest_new = (y['pop_age'][tipo][-1] / y['pop_age'][tipo][0]) * (1 - bq[-1]) * (1 + y['r']) * a_path[-1]

            death_dist = np.abs(bequest_new - bequest)

            bequest = (death_lam * bequest_new) + ((1 - death_lam) * bequest)

        # Now have the consistent initial asset holdings and the corresponding
        # path that goes with it using intratemporal and BC can calculate c,l

        c_path = np.zeros(life)
        l_path = np.zeros(life)

        for j in range(life - 1):
            y['bq'] = bq[j]
            if j <= work:
                c, l = agent_simul(a_path[j+1], a_path[j], y)
                c_path[j] = c
                l_path[j] = l
            else:
                y['version'] = 'retire'
                c, l = agent_simul(a_path[j+1], a_path[j], y)
                c_path[j] = c
                l_path[j] = l

        ######### Edit this section if you take retirement out ########

        y['version'] = 'retire'
        y['bq'] = bq[-1]
        c, l = agent_simul(a_grid[0], a_path[-1], y)

        c_path[-1] = c
        l_path[-1] = l

        return l_path, a_path, c_path, tipo

def gov_func(HH_res, args):

    taxes = args['taxes']

    total_revenue = 0

    gov_rev = {}
    gov_spe = {}

    for tipo in args['types_key']:

        l = HH_res[tipo][0]
        a = HH_res[tipo][1]
        c = HH_res[tipo][2]

        g_revenue = (args['wages'][tipo] * args['eff'][tipo] * taxes['labor tax'] * l) + (a * args['r'] * taxes['capital tax'])
        g_spending = taxes['transfers']

        g_revenue = np.sum(g_revenue)
        g_spending = np.sum(g_spending)

        gov_rev[tipo] = g_revenue
        gov_spe[tipo] = g_spending

        total_revenue += (gov_rev[tipo] - gov_spe[tipo])

    return total_revenue
