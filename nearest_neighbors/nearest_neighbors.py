import numpy as np
import sys
import matplotlib.pyplot as plt
from sim_code_do_not_change import *
from autocorrelation import *

# define parameters (k,a,N,m)
N = 9
m = np.full(N,12.)
k = np.zeros((N,N))
b = np.zeros((N,N))
for i in range(N):
    for j in range(i+1,N):
        k[i,j] = 303.1 # kcal/mol/A^2
        k[j,i] = k[i,j]
        b[i,j] = 2
        b[j,i] = -b[i,j]
        if j-i > 1:
            k[i,j] = 0
            k[j,i] = 0
print (k)
# define initial coords/vel (x,v)
x = np.zeros(N)
v = np.zeros(N)
for i in range(N):
    x[i] = 2.0*i

# define input to control sim (T,dt,Nstep,Nwrite)
T = 400.  # K
T *= 0.0019871571193690245 # kcal/mol 
gamma = 0.002 # fs^-1
gamma /= 0.02045482828087295 # MD type units
dt = 2.  # fs
dt *= 0.02045482828087295 # MD type units
c1 = np.exp(-gamma*dt)
c2 = np.sqrt(1-c1*c1)*np.sqrt(T)
Nstep = 3000001 
Nwrite = 3000 

allostery_simulation_1D(N,m,k,b,x,v,T,gamma,dt,c1,c2,Nstep,Nwrite)

#autocorrelation_particle1(N,dt)

