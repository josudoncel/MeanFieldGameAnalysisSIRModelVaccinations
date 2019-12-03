# -*- coding: utf-8 -*-
################################################################################################
#
# INPUT OF THE CODE: INITIAL PROPORTION OF SUSEPTIBLE AND INFECTED POPULATION (S0 AND I0, RESP.) 
#
# OUTPUT OF THE CODE: 
#     (a) evolution over time of the susceptible and infected population for the mean-field equilibrium
#     (b) 1: when at time zero the mean-field equilibrium vaccinates with maximum rate
#     (c) 0: when at time zero the mean-field equilibrium does not vaccinate
#
################################################################################################
   

import sys;
#import matplotlib.pyplot as plt;
import numpy as np;

c_V = 0.5;            # cost per unit time of vaccination
c_I = 36.5;           # cost per unit time of infection

beta = 73;            # infection rate
gamma = 36.5;         # recovery rate
vac_min=0;            # minimum vaccination rate
vac_max=10;           # maximum vaccination rate
T = 1;                           # time horizon
h = 0.001;
C = int(1.*T/h);                  # number of the intervals of the time grid
t  = np.linspace(0, T, C+1);      # time grid from 0 to T

my_policy = np.zeros(C+1);

# initial conditions
S0 = float(sys.argv[1]);              
I0 = float(sys.argv[2]);                 

j=0;
v=np.zeros(C+1);
x_S=np.zeros(C+1);
x_I=np.zeros(C+1);

for t_critical in xrange(C+1):
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

    v[j]= sum(1.*c_I*x_I*h+1.*c_V*x_S*my_policy*h);    
    j = j+1;

mymin=v[0];
mypos=0;
for i in xrange(1,C+1):
    if (v[i]<mymin):
        mymin = v[i];
        mypos = i;

print mypos;

for i in xrange(C):
        print x_S[i],x_I[i];
        if i<mypos:
            myvac=vac_max;
            my_policy[i]=vac_max;
        else:
            myvac=vac_min;
            my_policy[i]=vac_min;
        x_S[i+1] = x_S[i]+(-beta*x_I[i]*x_S[i] - myvac*x_S[i])*h;
        x_I[i+1] = x_I[i]+(beta*x_I[i]*x_S[i] - gamma*x_I[i])*h;    

if mypos>0:
    print 1;
else:
    print 0;

