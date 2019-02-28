import numpy as np
import matplotlib.pyplot as plt

Nstep = 2000001 
Nwrite= 2000


# plot
data = np.loadtxt("1D_allosterysim.dat")

pos = np.loadtxt("positions.dat")
nSteps = len(pos)
nSteps_range = range(nSteps)

t = np.array(range(int((Nstep-1)/Nwrite)))*2

plt.plot(t,pos[:,0],label='Particle 1')
plt.plot(t,pos[:,1],label='Particle 2')
plt.plot(t,pos[:,2],label='Particle 3')
plt.plot(t,pos[:,3],label='Particle 4')
plt.plot(t,pos[:,4],label='Particle 5')
plt.plot(t,pos[:,5],label='Particle 6')
plt.plot(t,pos[:,6],label='Particle 7')
plt.plot(t,pos[:,7],label='Particle 8')
plt.plot(t,pos[:,8],label='Particle 9')

plt.xlabel("Femtoseconds",size=14)
plt.ylabel("Angstroms",size=14)
plt.legend(loc='upper left')
plt.title("Positions vs. Time")
plt.grid()
plt.tight_layout()

plt.savefig("9N_positionvstime_w_centerofgeom.pdf",dpi=600,transparant=True)
plt.close()


# center of geometry motion
for ts in nSteps_range:
    center_of_geometry = np.mean(pos[ts])
    pos[ts] -= center_of_geometry
plt.hist(pos[:,0],histtype='barstacked',label='Particle 1',bins=500)
plt.hist(pos[:,1],histtype='step',label='Particle 2',bins=500)
plt.hist(pos[:,2],histtype='barstacked',label='Particle 3',bins=500)
plt.hist(pos[:,3],histtype='barstacked',label='Particle 4',bins=500)
plt.hist(pos[:,4],histtype='step',label='Particle 5',bins=500)
plt.hist(pos[:,5],histtype='barstacked',label='Particle 6',bins=500)
plt.hist(pos[:,6],histtype='barstacked',label='Particle 7',bins=500)
plt.hist(pos[:,7],histtype='step',label='Particle 8',bins=500)
plt.hist(pos[:,8],histtype='barstacked',label='Particle 9',bins=500)

plt.xlabel("Positions",size=14)
plt.ylabel("Counts",size=14)
plt.legend(loc='upper left')
plt.title("Position")
plt.grid()
plt.tight_layout()

plt.savefig("9N_positionhist.pdf",dpi=600,transparant=True)
plt.close()

plt.hist(data[:,0], histtype='step',label='Particle 1',bins=500)
plt.hist(data[:,3], histtype='step',label='Particle 2',bins=500)
plt.hist(data[:,6], histtype='step',label='Particle 3',bins=500)
plt.hist(data[:,9], histtype='step',label='Particle 4',bins=500)
plt.hist(data[:,12], histtype='step',label='Particle 5',bins=500)
plt.hist(data[:,15], histtype='step',label='Particle 6',bins=500)
plt.hist(data[:,18], histtype='step',label='Particle 7',bins=500)
plt.hist(data[:,21], histtype='step',label='Particle 8',bins=500)
plt.hist(data[:,24], histtype='step',label='Particle 9',bins=500)

plt.xlabel("Forces",size=14)
plt.ylabel("Counts",size=14)
plt.legend(loc='upper left')
plt.title("Forces")
plt.grid()
plt.tight_layout()

plt.savefig("9N_forceshist.pdf",dpi=600,transparant=True)
plt.close()

plt.hist(data[500:,2], histtype='step',label='Particle 1',bins=50)
plt.hist(data[500:,5], histtype='step',label='Particle 2',bins=50)
plt.hist(data[500:,8], histtype='step',label='Particle 3',bins=50)
plt.hist(data[500:,11], histtype='step',label='Particle 4',bins=50)
plt.hist(data[500:,14], histtype='step',label='Particle 5',bins=50)
plt.hist(data[500:,17], histtype='step',label='Particle 6',bins=50)
plt.hist(data[500:,20], histtype='step',label='Particle 7',bins=50)
plt.hist(data[500:,23], histtype='step',label='Particle 8',bins=50)
plt.hist(data[500:,26], histtype='step',label='Particle 9',bins=50)

plt.xlabel("Velocities",size=14)
plt.ylabel("Counts",size=14)
plt.legend(loc='upper left')
plt.title("Velocities")
plt.grid()
plt.tight_layout()

plt.savefig("9N_velocitieshist.pdf",dpi=600,transparant=True)
plt.close()


t = np.array(range(int((Nstep-1)/Nwrite)))*2

plt.plot(t,pos[:,0],label='Particle 1')
plt.plot(t,pos[:,1],label='Particle 2')
plt.plot(t,pos[:,2],label='Particle 3')
plt.plot(t,pos[:,3],label='Particle 4')
plt.plot(t,pos[:,4],label='Particle 5')
plt.plot(t,pos[:,5],label='Particle 6')
plt.plot(t,pos[:,6],label='Particle 7')
plt.plot(t,pos[:,7],label='Particle 8')
plt.plot(t,pos[:,8],label='Particle 9')

plt.ylim(-10,10)
plt.xlabel("Femtoseconds",size=14)
plt.ylabel("Angstroms",size=14)
plt.legend(loc='upper left')
plt.title("Positions vs. Time")
plt.grid()
plt.tight_layout()

plt.savefig("9N_positionvstime.pdf",dpi=600,transparant=True)
plt.close()

t = np.array(range(int((Nstep-1)/Nwrite)-500))*2

plt.plot(t,data[500:,2],label='Particle 1')
plt.plot(t,data[500:,5],label='Particle 2')
plt.plot(t,data[500:,8],label='Particle 3')
plt.plot(t,data[500:,11],label='Particle 4')
plt.plot(t,data[500:,14],label='Particle 5')
plt.plot(t,data[500:,17],label='Particle 6')
plt.plot(t,data[500:,20],label='Particle 7')
plt.plot(t,data[500:,23],label='Particle 8')
plt.plot(t,data[500:,26],label='Particle 9')

plt.xlabel("Femtoseconds",size=14)
plt.ylabel("Velocity",size=14)
plt.legend(loc='upper left')
plt.title("Velocity vs. Time")
plt.grid()
plt.tight_layout()

plt.savefig("9N_velocitiesvstime.pdf",dpi=600,transparant=True)
plt.close()

