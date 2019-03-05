import numpy as np
import sys
import matplotlib.pyplot as plt
from fast_sim_code import *


# define input to control sim (T,dt,Nstep,Nwrite)
T = 303.  # K
gamma = 0.002 # fs^-1
dt = 2.  # fs
Nstep = 201
Nwrite = 2

# read initial coordinates
pos = np.loadtxt("pos_initial.dat")
m = 12.
N = pos.shape[0]
# read force field parameters as NxN matrices
k = -np.loadtxt("IGPS_hessian.dat")
b = np.loadtxt("IGPS_pair_dist.dat")
# sanity check
if N != k.shape[0] or N != b.shape[0]:
    print("size of position array does not agree with force field matrices")
    sys.exit(0)
# debug using two particles
#N = 2
#pos = pos[:2,:]
#k = k[:2,:2]
#b = b[:2,:2]
# declare force and velocity arrays
f = np.zeros((N,3))
v = np.zeros((N,3))

# convert Langevin thermostat stuff
T *= 0.0019871571193690245 # kcal/mol 
gamma /= 0.02045482828087295 # MD type units
dt *= 0.02045482828087295 # MD type units
c1 = np.exp(-gamma*dt)
c2 = np.sqrt((1-c1**2)*T/m)

# make array of nodes and edges (atoms and force constants)
edges = []
nEdges = 0

# run simulation
allostery_simulation_3D(N,m,k,b,pos,v,T,f,gamma,dt,c1,c2,Nstep,Nwrite,edges,nEdges)

