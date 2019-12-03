# -*- coding: utf-8 -*-

################################################################################################
#
# INPUT OF THE CODE: INITIAL PROPORTION OF SUSEPTIBLE AND INFECTED POPULATION (S0 AND I0, RESP.) 
#
# OUTPUTS OF THE CODE: 
#     (a) evolution over time of the susceptible and infected population for the mean-field equilibrium
#     (b) 1: when at time zero the mean-field equilibrium vaccinates with maximum rate
#     (c) 0: when at time zero the mean-field equilibrium does not vaccinate
#
################################################################################################
   

import sys;
import numpy as np;
#import matplotlib.pyplot as plt;

plt.close();        

def equilibrium(c_V):
    h=1.*T/C;
    if (1./h < max(beta,gamma,vac_max)):
        return -1;
    t_critical=0;#randint(0,C); 
    thr=-1;
    res_a=np.zeros(C+1);

    while (t_critical != thr):
        
        x_S = np.zeros(C+1); x_I = np.zeros(C+1);
        my_policy = np.zeros(C+1);
        x_S[0]=S0; x_I[0]=I0;
        for i in xrange(C):
            if i<t_critical:
                myvac=vac_max;
                my_policy[i]=vac_max;
            else:
                myvac=vac_min;
                my_policy[i]=vac_min;
            x_S[i+1] = x_S[i]+(-beta*x_I[i]*x_S[i] - myvac*x_S[i])*h;
            x_I[i+1] = x_I[i]+(beta*x_I[i]*x_S[i] - gamma*x_I[i])*h;    

        J_s = np.zeros(C+1);
        J_i = np.zeros(C+1);
        J_r = np.zeros(C+1);
        J_v = np.zeros(C+1);

        for count in xrange(C):
            t_ = C - count ;
            slope = c_V - J_s[t_];
            if slope < 0:
                res_a[t_]=vac_max;
                p_s_s=1 - 1.*beta*x_I[t_-1]*h - 1.*vac_max*h;
                p_s_i=1.*beta*x_I[t_-1]*h; 
                J_s[t_-1]= vac_max * c_V *h + p_s_s * J_s[t_] + p_s_i * J_i[t_];
            else:
                res_a[t_]=vac_min;
                p_s_s=1 - beta*x_I[t_-1]*h - vac_min*h;
                p_s_i=beta*x_I[t_-1]*h; 
                J_s[t_-1]= vac_min * c_V *h + p_s_s * J_s[t_] + p_s_i * J_i[t_];   
            p_i_s=0; 
            p_i_i=1-gamma*h; 
            J_i[t_-1] = c_I*h + p_i_s * J_s[t_] + p_i_i * J_i[t_];   
            J_r[t_-1]=J_r[t_];  J_v[t_-1]=J_v[t_]; 

        thr=0;
        for i in xrange(C+1):
            if J_s[i]<c_V:
                thr=i;
                break;
        
        if t_critical > thr:
            t_critical = t_critical - 1;
        if t_critical < thr:
            t_critical = t_critical + 1;  
            
        #print t_critical, thr, J_s;


    return thr;

c_V = 0.5;           # cost per unit time of vaccination
c_I = 36.5;               # cost per unit time of infection

beta = 73;            # infection rate
gamma = 36.5;         # recovery rate
vac_min=0;            # minimum vaccination rate
vac_max=10;           # maximum vaccination rate
T = 1;                           # time horizon
h = 0.001;
C = int(1.*T/h);                  # number of the intervals of the time grid
t  = np.linspace(0, T, C+1);      # time grid from 0 to T


                         
# initial conditions
S0 = float(sys.argv[1]);            
I0 = float(sys.argv[2]);                 

x_S = np.zeros(C+1);
x_I = np.zeros(C+1);
t_critical = equilibrium(c_V);

my_policy = np.zeros(C+1);
x_S[0]=S0; x_I[0]=I0;
for i in xrange(C):
    print x_S[i],x_I[i]
    if i<t_critical:
        myvac=vac_max;
        my_policy[i]=vac_max;
    else:
        myvac=vac_min;
        my_policy[i]=vac_min;
    x_S[i+1] = x_S[i]+(-beta*x_I[i]*x_S[i] - myvac*x_S[i])*h;
    x_I[i+1] = x_I[i]+(beta*x_I[i]*x_S[i] - gamma*x_I[i])*h; 
    
if t_critical>0:
    print 1;
else:
    print 0;



