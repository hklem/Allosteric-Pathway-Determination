import numpy as np
import sys
import matplotlib.pyplot as plt
from fast_sim_code import *

# define parameters (k,b,N,m)
N = 454
#m = np.full(N,12.)
m = 12.
k = np.loadtxt("IGPS_hessian.dat")
D = len(k)
b = np.loadtxt("IGPS_pair_dist.dat")


# define initial coords/vel (xyz,v)
pos = np.loadtxt("pos_initial.dat")
x = np.zeros(N)
y = np.zeros(N)
z = np.zeros(N)
r = np.zeros((N,3))
# make array of nodes and edges (atoms and force constants)
edges = []
nEdges = 0

for i in range(N):
    x[i] = pos[i,0]
    y[i] = pos[i,1]
    z[i] = pos[i,2]
    r[i] = pos[i,:]
v = np.zeros((N,3))

# define input to control sim (T,dt,Nstep,Nwrite)
T = 303.  # K
T *= 0.0019871571193690245 # kcal/mol 
gamma = 0.002 # fs^-1
gamma /= 0.02045482828087295 # MD type units
dt = 2.  # fs
dt *= 0.02045482828087295 # MD type units
c1 = np.exp(-gamma*dt)
#c2 = np.sqrt(1-c1*c1)*np.sqrt(T)
c2 = np.sqrt((1-c1**2)*T/m)
Nstep = 21 
Nwrite = 2
f = np.zeros((N,3))

allostery_simulation_3D(N,m,k,b,r,v,T,f,gamma,dt,c1,c2,Nstep,Nwrite,edges,nEdges)

