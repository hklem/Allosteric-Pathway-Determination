# PREAMBLE

import numpy as np
import matplotlib.pyplot as plt
from pearson_correlation import *

def allostery_simulation_1D(N,m,k,b,x,v,T,gamma,dt,c1,c2,Nstep,Nwrite):
    # loop over timestep (calc forces and propagate)
    outputfile = open("1D_allosterysim.dat","w")
    positions = open("positions.dat","w")

    var_string = "#"

    for it in range(Nstep):
        print(it)
        f = np.zeros(N)
        var_string = ""
        pos_string = ""
        for i in range(N):
            # calc forces
            for j in range(i+1,N):
                f_ij = -k[i,j]*(x[i]-x[j]-b[i,j])
                f[i] += f_ij
                f[j] -= f_ij

        # langevin
        for i in range(N):
            v[i] += (f[i]/m[i])*dt/2. #B
            x[i] += v[i]*dt/2. #A
            v[i] = c1*v[i] + np.random.normal()*(c2/np.sqrt(m[i])) #O
            x[i] += v[i]*dt/2.     
            
        f = np.zeros(N)
        for i in range(N):
            # calc forces
            for j in range(i+1,N):
                f_ij = -k[i,j]*(x[i]-x[j]-b[i,j])
                f[i] += f_ij
                f[j] -= f_ij
            v[i] += (f[i]/m[i])*dt/2. #B
            var_string += str(f[i])+"  "+str(x[i])+"  "+str(v[i])+"  "
            pos_string += str(x[i])+"  "    
        if it%Nwrite == 0:
            if it > 0:
                outputfile.write(var_string+"\n")
                positions.write(pos_string+"\n")

    outputfile.close()
    positions.close()

    print(var_string)
    print(pos_string)
    
    # pearson_correlation calculation
    trajectory_data = np.loadtxt("positions.dat")
    node_correlation, node_variance, node_average = pearson_correlation(trajectory_data)

    c = np.zeros((N,N)) 
    for i in range(N):
        for j in range(i,N):
            c[i,j] = node_correlation[i,j]*np.sqrt(node_variance[i]*node_variance[j])
            c[j,i] = c[i,j]

    np.savetxt("covariance_nearest_neighbors.dat",c)

