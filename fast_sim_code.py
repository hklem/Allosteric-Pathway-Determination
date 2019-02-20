import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as MDA
import numba
from numba import jit

@jit
def iterations(N,m,k,b,r,v,T,f,gamma,dt,c1,c2,Nstep,Nwrite,edges,nEdges):

    for i in range(N):
        # thermostat
        v[i] += (f[i]/m)*dt/2. #B
        r[i] += v[i]*dt/2. #A
        v[i] = c1*v[i] + np.random.normal(size=3)*c2 #O
        r[i] += v[i]*dt/2.     
        
    f = np.zeros((N,3))
        
    for edge in edges:
        i = edge[0]
        j = edge[1]
        k_ij = edge[2]
        b_ij = edge[3]
        # recalc forces
        f_ij = -k_ij*(1-(b_ij/(np.linalg.norm(r[i]-r[j])))*(r[j]-r[i]))
        f[i] -= f_ij
        f[j] += f_ij
    for i in range(N):
        v[i] += (f[i]/m)*dt/2. #B
@jit
def initial_forces(N,r,f,edges):
    f = np.zeros((N,3))
    for edge in edges:
        i = edge[0]
        j = edge[1]
        k_ij = edge[2]
        b_ij = edge[3]
        # calc forces
        f_ij = -k_ij*(1-(b_ij/(np.linalg.norm(r[i]-r[j])))*(r[j]-r[i]))
        f[i] -= f_ij
        f[j] += f_ij

def allostery_simulation_3D(N,m,k,b,r,v,T,f,gamma,dt,c1,c2,Nstep,Nwrite,edges,nEdges):
    for i in range(N-1):
        for j in range (i+1,N):
            if -k[i,j] > 1e-5:
                #edges.append([])
                #edges[nEdges].append(i)
                #edges[nEdges].append(j)
                #edges[nEdges].append(k[i,j])
                #edges[nEdges].append(b[i,j])
                edges.append([i,j,k[i,j],b[i,j]])
                nEdges += 1
            
  # loop over timestep (calc forces and propagate)
    forces_xyz = open("forces_xyz.dat","w")
    velocities_xyz = open("velocities_xyz.dat","w")
   
    u = MDA.Universe("alignedCasLastFramed.pdb")
    allu = u.select_atoms("all")
    positions_xyz = MDA.Writer("positions.dcd",N)
    velocities_xyz = MDA.Writer("velocities.dcd",N)
    forces_xyz = MDA.Writer("forces.dcd",N)
  
   
   # define initial coords/vel (xyz,v)
    initial_forces(N,r,f,edges)
    for it in range(Nstep):
        print(it)
        iterations(N,m,k,b,r,v,T,f,gamma,dt,c1,c2,Nstep,Nwrite,edges,nEdges)
        if it%Nwrite == 0:
            if it > 0:
                allu.positions = r        
                positions_xyz.write(allu)
                allu.positions = f
                forces_xyz.write(allu)
                allu.positions = v
                velocities_xyz.write(allu)
    positions_xyz.close()
    forces_xyz.close()
    velocities_xyz.close()
 
 
 


