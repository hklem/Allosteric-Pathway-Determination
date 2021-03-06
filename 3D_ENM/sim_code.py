import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as MDA
import numba
from numba import jit

@jit
def iterations(N,m,k,b,r,v,T,f,gamma,dt,c1,c2,Nstep,Nwrite):
    for i in range(N):
        # thermostat
        v[i] += (f[i]/m)*dt/2. #B
        r[i] += v[i]*dt/2. #A
        v[i] = c1*v[i] + np.random.normal(size=3)*c2 #O
        r[i] += v[i]*dt/2.     
        
    # recalc forces
    f = np.zeros((N,3))
    for i in range(N):
        for j in range(i+1,N):
            if k[i,j] > 1e-3:
                f_ij = -k[i,j]*(1-(b[i,j]/(np.linalg.norm(r[i]-r[j]))))*(r[j]-r[i])
                f[i] -= f_ij
                f[j] += f_ij
    for i in range(N):
        v[i] += (f[i]/m)*dt/2. #B
    return r, v, f
@jit
def initial_forces(N,r,f,k,b):
    f = np.zeros((N,3))
    for i in range(N):
        for j in range(i+1,N):
            if k[i,j] > 1e-3:
                # calc forces
                f_ij = -k[i,j]*(1-(b[i,j]/(np.linalg.norm(r[i]-r[j]))))*(r[j]-r[i])
                f[i] -= f_ij
                f[j] += f_ij
    return f

def allostery_simulation_3D(N,m,k,b,r,v,T,f,gamma,dt,c1,c2,Nstep,Nwrite):
  # loop over timestep (calc forces and propagate)
    u = MDA.Universe("alignedCasLastFramed.pdb")
    allu = u.select_atoms("all")
    positions_xyz = MDA.Writer("positions.dcd",N)
    velocities_xyz = MDA.Writer("velocities.dcd",N)
    forces_xyz = MDA.Writer("forces.dcd",N)
   
   # define initial coords/vel (xyz,v)
    f = initial_forces(N,r,f,k,b)
    for it in range(Nstep):
        print(it)
        r, v, f = iterations(N,m,k,b,r,v,T,f,gamma,dt,c1,c2,Nstep,Nwrite)
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
 
