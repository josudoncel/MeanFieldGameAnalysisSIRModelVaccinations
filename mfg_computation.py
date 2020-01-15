"""
This module contains all the code needed to compute the mean field equilibrium and the social optimum
"""
import numpy as np


class Parameters:
    """
    Contains all parameters of a model
    """

    def __init__(self):
        """
        Default values (used in the paper)
        """
        self.gamma = 73.0     # infection rate
        self.rho = 36.5    # recovery rate
        self.vac_min = 0.0   # minimum vaccination rate
        self.theta = 10.0  # maximum vaccination rate
        self.T = 1           # time horizon

        self.c_V = 0.5       # cost per unit time of vaccination
        self.c_I = 36.5      # cost per unit time of infection

def compute_switching_curve(parameters, pol):
    M = len(pol)
    theta = parameters.theta; vac_min = parameters.vac_min 
    sw_x = []; sw_y = []; 
    for ms in range(M):
        for mi in range(M):
            if pol[ms,mi] == theta and pol[ms,mi-1] == vac_min:
                sw_x.append(ms/M); sw_y.append(mi/M)
    return sw_x, sw_y
        
        
def mf_equilibrium(parameters, S0, I0, C=100000):
    """
    We compute the mean-field equilibrium as a fixed-point problem. The algorithm stops if
    after n_loops_max the fixed point has not been found.
    """
    time_horizon = int(parameters.T*C)
    gamma, rho, vac_min, theta = parameters.gamma, parameters.rho, parameters.vac_min, parameters.theta
    c_V, c_I = parameters.c_V, parameters.c_I
    best_pol = np.zeros(time_horizon+1)
    value_s = np.zeros(time_horizon+1)
    value_i = np.zeros(time_horizon+1)

    def best_response_threshold(thr_mass):
        x_S = np.zeros(time_horizon+1)
        x_I = np.zeros(time_horizon+1)
        policy = np.zeros(time_horizon+1)
        x_S[0] = S0
        x_I[0] = I0
        for i in range(time_horizon):
            if i < thr_mass:
                policy[i] = theta
            else:
                policy[i] = vac_min
            x_S[i+1] = x_S[i]+(-gamma*x_I[i]*x_S[i] - policy[i]*x_S[i])/C
            x_I[i+1] = x_I[i]+(gamma*x_S[i] - rho)*x_I[i]/C

        for t_ in range(time_horizon, 0, -1):
            if c_V - value_s[t_] < 0:
                best_pol[t_-1] = theta
                p_s_s = 1 - 1.*gamma*x_I[t_-1]/C - 1.*theta/C
                p_s_i = 1.*gamma*x_I[t_-1]/C
                value_s[t_-1] = theta * c_V / C + p_s_s * \
                    value_s[t_] + p_s_i * value_i[t_]
            else:
                best_pol[t_-1] = vac_min
                p_s_s = 1 - 1.*gamma*x_I[t_-1]/C - vac_min/C
                p_s_i = gamma*x_I[t_-1]/C
                value_s[t_-1] = vac_min * c_V/C + p_s_s * \
                    value_s[t_] + p_s_i * value_i[t_]
            p_i_i = 1-rho/C
            value_i[t_-1] = c_I/C + p_i_i * value_i[t_]

        for i in range(time_horizon+1):
            if value_s[i] < c_V:
                thr_pl0 = i
                break
        return thr_pl0

    # Now we first find a threshold such that the BR is smaller than the threshold
    thr_mass = 0
    while best_response_threshold(thr_mass) > thr_mass:
        thr_mass = max(1, thr_mass*2)

    # Then we apply a dicotomic search
    thr_mass = int(thr_mass/2)
    delta_threshold = int(np.ceil(thr_mass/2))
    for i in range(20):
        br_threshold = best_response_threshold(thr_mass)
        if br_threshold > thr_mass:
            thr_mass += delta_threshold
        elif br_threshold < thr_mass:
            thr_mass -= delta_threshold
        delta_threshold = min(thr_mass, int(np.ceil(delta_threshold/2)))
        if delta_threshold <= 1:
            break

    cost_mfe = value_s[0]*S0+value_i[0]*I0
    return cost_mfe, br_threshold/C


def mf_optimum(parameters, S0, I0, C=100000):
    """
    This method uses a divide and conquer approach to compute the minimum value.
    It assumes that the optimal policy is a bang-bang policy with a unique jump time"""
    if S0 == 0:
        t = 0.0
        return S0, t
    time_horizon = int(parameters.T*C)
    gamma, rho, vac_min, theta = parameters.gamma, parameters.rho, parameters.vac_min, parameters.theta
    c_V, c_I = parameters.c_V, parameters.c_I

    my_policy = np.zeros(time_horizon+1)
    v = np.zeros(time_horizon+1)
    x_S = np.zeros(time_horizon+1)
    x_I = np.zeros(time_horizon+1)

    def compute_cost(t_critical):
        x_S[0] = S0
        x_I[0] = I0
        if v[t_critical] > 0:
            return
        for i in range(time_horizon):
            my_policy[i] = theta if i < t_critical else vac_min
            x_S[i+1] = x_S[i]+(-gamma*x_I[i]*x_S[i] - my_policy[i]*x_S[i])/C
            x_I[i+1] = x_I[i]+(gamma*x_I[i]*x_S[i] - rho*x_I[i])/C

        v[t_critical] = sum(1.*c_I*x_I/C+1.*c_V*x_S*my_policy/C)

    t_critical = int(C/2)
    deltaT = int(C/4)
    for i in range(20):
        for t in [t_critical, t_critical+1]:
            compute_cost(t)
        if v[t_critical+1] > v[t_critical]:
            t_critical -= deltaT
        else:
            t_critical += deltaT
        deltaT = min(t_critical, int(np.ceil(deltaT/2)))
        if deltaT == 0:
            break

    cost_opt = v[t_critical]
    return cost_opt, t_critical/C


def trans_prob(m_s, m_i, m_pi, pl0_pi, N, C, PAR):
    gamma = PAR.gamma; rho = PAR.rho
    p_x = np.zeros(3)
    p_m = np.zeros(4)
    # transition probabilities of Player 0
    p_x[0] = 1. * gamma / C * m_i / N # proba of infection
    p_x[1] = 1. * rho / C             # proba of recovering 
    p_x[2] = 1. * pl0_pi / C          # proba of vaccination
    # transition probabilities of the mass
    p_m[0] = 1. * gamma * N / C * m_s / N * m_i / N # proba of infection when Player0=Susceptible
    p_m[1] = 1. * m_i / N * rho * N / C             # proba of recovering
    p_m[2] = 1. * m_s / N * m_pi * N / C            # proba of vaccination 
    p_m[3] = 1. * gamma * N / C * m_s / N * (m_i+1) / N # proba of infection when Player0=Infected
    return p_x, p_m

def uniformization_constant(N):
    return T*N*100 #gamma / N + rho + theta + N * (gamma/4 + rho + theta)* N 

def best_response_n_players(m_pol, C, N, params):

    theta = params.theta
    gamma = params.gamma
    rho = params.rho
    c_I = params.c_I
    c_V = params.c_V

    if C < N * max(theta, gamma, rho): raise Exception('Higher value of C required')

    # Pl0 cost at susceptible and infected state: Js and Ji
    Js = np.zeros((C + 1, N + 1, N + 1))
    Ji = np.zeros((C + 1, N + 1, N + 1))
    # Best response of Player 0
    BR_pol = np.zeros((C + 1, N + 1, N + 1))

    for t in range(C, 0, -1):
        for ms in range(0, N + 1):
            for mi in range(0, N - ms + 1):
                if c_V - Js[t, ms, mi] < 0:
                    BR_pol[t - 1, ms, mi] = theta
                else:
                    BR_pol[t - 1, ms, mi] = 0
                px, pm = trans_prob(ms, mi, m_pol[t - 1, ms, mi], BR_pol[t - 1, ms, mi], N, C, params)
                Js[t - 1, ms, mi] = (
                    c_V * BR_pol[t - 1, ms, mi] / C +                               # pl0 instantaneous cost at susceptible state
                    (px[0] * Ji[t, ms, mi] if mi > 0 else 0) +                      # pl0 infection
                    (px[2] * 0) +                                                   # pl0 vaccination: the cost when it is recovered is 0
                    (pm[0] * Js[t, ms - 1, mi + 1] if ms > 0 else 0) +              # infection of mass (when Player0=S)
                    (pm[1] * Js[t, ms, mi - 1] if mi > 0 else 0) +                  # recovery of mass
                    (pm[2] * Js[t, ms - 1, mi] if ms > 0 else 0) +                  # vaccination of mass
                    (1 - px[0] - px[2] - pm[0] - pm[1] - pm[2]) * Js[t, ms, mi])    # no transition
                #print((1 - px[0] - px[2] - pm[0] - pm[1] - pm[2]))
                Ji[t - 1, ms, mi] = (
                    c_I / C +                                               # pl0 instantaneous cost at infected state
                    (px[1] * 0) +                                           # pl0 recovery: the cost when it is recovered is 0
                    (pm[3] * Ji[t, ms - 1, mi + 1] if ms > 0 else 0) +      # infection of mass (when Player0=I)
                    (pm[1] * Ji[t, ms, mi - 1] if mi > 0 else 0) +          # recovery of mass
                    (pm[2] * Ji[t, ms - 1, mi] if ms > 0 else 0) +          # vaccination of mass
                    (1 - px[1] - pm[3] - pm[1] - pm[2]) * Ji[t, ms, mi])    # no transition
    return BR_pol, Js, Ji



def fileName(N):
    return 'data/equi_N{}.npy'.format(N)

def initialization(params,N):
    try:
        equilibrium_N30, Js, Ji = np.load(fileName(30))
        C = params.T*N*100 #uniformization_constant(N)
        m_pol = np.zeros((C+1,N+1,N+1))
        for t in range(C+1):
            for ms in range(N+1):
                for mi in range(N+1):
                    m_pol[t,ms,mi] = equilibrium_N30[ int(t/N*30), int(ms/N*30), int(mi/N*30)]
        return m_pol
    except Exception as e:
        print(e)
        m_pol = np.ones((C + 1, N + 1, N + 1))*0
        for ms in range(0,N+1):
            for mi in range(0,N+1-ms):
                if 2*ms+mi >= 0.7*N:
                    m_pol[:,ms,mi] = theta
        return m_pol

def get_equilibrium_n_players(N):
    try:
        BR_pol, Js, Ji = np.load(fileName(N))
        return BR_pol, Js, Ji
    except Exception as e:
        print(e)
    
def compute_equilibrium_n_players(parameters, N):

    T = parameters.T
    C = T*N*100 #uniformization_constant(N)
    # the policy of the mass is never vaccinate
    m_pol = initialization(parameters,N)
    
    # the policy of the mass is to always vaccinate with max rate
    # m_pol=theta*np.ones((C+1,N+1,N+1));

    BR_pol, Js, Ji = best_response_n_players(m_pol, C, N, parameters)
    # print(BR_pol[0])

    iteration_numbers = 0
    while True:
        number_of_differences = np.sum(1-np.isclose(BR_pol, m_pol))
        if number_of_differences == 0:
            print("Equilibrium found. Equilibrium policy at t=0:")
            print(BR_pol[0])
            break
        else:
            print("Ite{:3}: BR_pol and m_pol do not coindice for {:5} values (diff={:7.2f}) Updating m_pol...".format(
                iteration_numbers, number_of_differences, np.sum(BR_pol - m_pol)))
            if np.random.rand() >= 0.1:
                m_pol = BR_pol
            else:
                print('*',end='')
                m_pol = BR_pol*0.2+m_pol*0.8 # This code is to try to 
            BR_pol, Js, Ji = best_response_n_players(m_pol, C, N, parameters)
        iteration_numbers += 1
        
    np.save(fileName(N), np.array([BR_pol, Js, Ji]))
    return BR_pol, Js, Ji

def get_globaloptimum_n_players(N):
    try:
        BR_pol, J = np.load('data/globalopt_N{}.npy'.format(N))
        return BR_pol, J
    except Exception as e:
        print(e)

def compute_globaloptimum_n_players(parameters,N):

    T = parameters.T
    theta = parameters.theta
    gamma = parameters.gamma
    rho = parameters.rho
    c_I = parameters.c_I
    c_V = parameters.c_V
    C = T*N*100 #uniformization_constant(N)
    if C < N * max(theta, gamma, rho): raise Exception('Higher value of C required')
    J = np.zeros((N+1,N+1,C+1))
    bestpol = np.zeros((N+1,N+1,C+1))
    for t in range(C,0,-1):
        # CASE: no susceptibles and no infected. The bestpol is never vaccinate 
        J[0,0,t-1] = 0
        bestpol[0,0,t-1] = 0
        # CASE: no susceptibles, but there are infected. The bestpol is undefined (we set 'never vaccinate')
        for Mi in range(1,N+1):
            Ms = 0
            J[Ms,Mi,t-1]= c_I * Mi / N / C \
                        + rho / C * Mi * J[Ms,Mi-1,t] \
                        + (1 - rho / C * Mi) * J[Ms,Mi,t]
            bestpol[Ms,Mi,t-1] = 0
        for Ms in range(1,N+1):
            # CASE: susceptibles, but no infected. The bestpol is never vaccinate
            J[Ms,0,t-1] = 0
            bestpol[Ms,Mi,t-1] = 0
            # CASE: susceptibles and infected.
            for Mi in range(1,N-Ms+1):
                if c_V/N +  J[Ms-1,Mi,t] -  J[Ms,Mi,t] < 0:
                    # bestpol is to vaccinate with maximum rate
                    J[Ms,Mi,t-1] = (c_V * theta * Ms / N + c_I * Mi / N) / C \
                                + theta / C * Ms  * J[Ms-1,Mi,t] \
                                + rho / C * Mi * J[Ms,Mi-1,t]\
                                + gamma / C * Mi * Ms / N * J[Ms-1,Mi+1,t]\
                                +(1-theta/C*Ms-rho/C*Mi-gamma/C*Mi*Ms/N)*J[Ms,Mi,t];
                    bestpol[Ms,Mi,t-1] = theta
                else:
                    # bestpol is not to vaccinate
                    J[Ms,Mi,t-1]=c_I * Mi / N / C \
                                + rho / C * Mi * J[Ms,Mi-1,t]\
                                + gamma / C * Mi * Ms / N *J[Ms-1,Mi+1,t]\
                                + (1 - rho / C * Mi - gamma / C * Mi * Ms / N)*J[Ms,Mi,t];
                    #print((1-rho/C*Mi/N-gamma/C*Mi/N*Ms/N))
                    bestpol[Ms,Mi,t-1] = 0
                    
    np.save('data/globalopt_N{}_players.npy'.format(N), np.array([bestpol, J]))
    return bestpol, J

# To test the code, uncomment the following lines
#PARAMS = Parameters()
#print('cost of mfe=', mf_equilibrium(PARAMS, 0.4, 0.4, 100000))
#print('cost of opt=', mf_optimum(PARAMS, 0.4, 0.4, 100000))


