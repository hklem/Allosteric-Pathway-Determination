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

b[0,2] = 4
b[2,0] = -4
b[1,3] = 4
b[3,1] = -4
b[0,3] = 6
b[3,0] = -6
b[5,7] = 4
b[7,5] = -4
b[5,8] = 6
b[8,5] = -6
b[6,8] = 4
b[8,6] = -4

k[0,4] = 0
k[4,0] = 0
k[1,4] = 0
k[4,1] = 0
k[2,4] = 0
k[4,2] = 0
k[4,6] = 0
k[6,4] = 0
k[4,7] = 0
k[7,4] = 0
k[4,8] = 0
k[8,4] = 0

k[0,5] = 0
k[5,0] = 0
k[0,6] = 0
k[6,0] = 0
k[0,7] = 0
k[7,0] = 0
k[0,8] = 0
k[8,0] = 0
k[1,5] = 0
k[5,1] = 0
k[1,6] = 0
k[6,1] = 0
k[1,7] = 0
k[7,1] = 0
k[1,8] = 0
k[8,1] = 0

k[2,5] = 0
k[5,2] = 0
k[2,6] = 0
k[6,2] = 0
k[2,7] = 0
k[7,2] = 0
k[2,8] = 0
k[8,2] = 0


k[3,5] = 0
k[5,3] = 0
k[3,6] = 0
k[6,3] = 0
k[3,7] = 0
k[7,3] = 0
k[3,8] = 0
k[8,3] = 0

print(k)
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
Nstep = 2000001 
Nwrite = 2000 

allostery_simulation_1D(N,m,k,b,x,v,T,gamma,dt,c1,c2,Nstep,Nwrite)

autocorrelation_particle1(N,dt)

