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
        
    for edge in range(nEdges):
        i = int(edges[edge,0])
        j = int(edges[edge,1])
        k_ij = edges[edge,2]
        b_ij = edges[edge,3]
        # recalc forces
        f_ij = -k_ij*(1-(b_ij/(np.linalg.norm(r[i]-r[j]))))*(r[j]-r[i])
        f[i] -= f_ij
        f[j] += f_ij
    for i in range(N):
        v[i] += (f[i]/m)*dt/2. #B
    return r, v, f

@jit
def initial_forces(N,r,edges):
    f = np.zeros((N,3))
    for edge in range(edges.shape[0]):
        i = int(edges[edge,0])
        j = int(edges[edge,1])
        k_ij = edges[edge,2]
        b_ij = edges[edge,3]
        # calc forces
        f_ij = -k_ij*(1-(b_ij/(np.linalg.norm(r[i]-r[j]))))*(r[j]-r[i])
        f[i] -= f_ij
        f[j] += f_ij
    return f

def allostery_simulation_3D(N,m,k,b,r,v,T,f,gamma,dt,c1,c2,Nstep,Nwrite,edges,nEdges):
    for i in range(N-1):
        for j in range (i+1,N):
            if k[i,j] > 1e-3:
                #edges.append([])
                #edges[nEdges].append(i)
                #edges[nEdges].append(j)
                #edges[nEdges].append(k[i,j])
                #edges[nEdges].append(b[i,j])
                edges.append([i,j,k[i,j],b[i,j]])
                nEdges += 1
    edges = np.asarray(edges)
  # loop over timestep (calc forces and propagate)
    f_xyz = open("sim_data/forces_xyz.dat","w")
    v_xyz = open("sim_data/velocities_xyz.dat","w")
    r_xyz = open("sim_data/positions_xyz.dat","w")
   
    u = MDA.Universe("alignedCasLastFramed.pdb")
#    u = MDA.Universe("2particles.pdb")
    allu = u.select_atoms("all")
    positions_xyz = MDA.Writer("sim_data/positions.dcd",N)
    velocities_xyz = MDA.Writer("sim_data/velocities.dcd",N)
    forces_xyz = MDA.Writer("sim_data/forces.dcd",N)
  
   
   # define initial coords/vel (xyz,v)
    f = initial_forces(N,r,edges)
    for it in range(Nstep):
        print(it)
        r, v, f = iterations(N,m,k,b,r,v,T,f,gamma,dt,c1,c2,Nstep,Nwrite,edges,nEdges)
        if it%Nwrite == 0:
            if it > 0:
                allu.positions = r        
                positions_xyz.write(allu)
                allu.positions = f
                forces_xyz.write(allu)
                allu.positions = v
                velocities_xyz.write(allu)
                for atom in range(N):
                    f_xyz.write("%16.8f %16.8f %16.8f\n" % (f[atom,0],f[atom,1],f[atom,2]))
                    v_xyz.write("%16.8f %16.8f %16.8f\n" % (v[atom,0],v[atom,1],v[atom,2]))
                    r_xyz.write("%16.8f %16.8f %16.8f\n" % (r[atom,0],r[atom,1],r[atom,2]))
    positions_xyz.close()
    forces_xyz.close()
    velocities_xyz.close()
    r_xyz.close()
    f_xyz.close()
    v_xyz.close()
 
 
 


