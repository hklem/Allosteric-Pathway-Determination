import numpy as np
import sys
import matplotlib.pyplot as plt

# define parameters (k,a,N,m)
N = 9
m =  np.full(N,12.)
k = np.zeros((N,N))
b = np.zeros((N,N))
for i in range(N):
        for j in range(i+1,N): 
                k[i,j] = 303.1
                k[j,i] = k[i,j]
                b[i,j] = 2.
                b[j,i] = -b[i,j]
print("parameters good")

# define initial coords/vel (x,v)
x = np.zeros(N)
v = np.zeros(N)
for i in range(N):
        x[i] = 1.5*i
print("initial coords/vel good")

# define input to control sim (T,dt,Nstep,Nwrite)
T = 400.  # K
T *= 0.0019871571193690245 # kcal/mol
gamma = 1
dt = 2.  # fs
dt *= 0.02045482828087295 # MD type units
c1 = np.exp(-gamma*dt)
c2 = np.sqrt(1-c1*c1)*np.sqrt(T)
Nstep = 10000
Nwrite = 100
print ("input good")

# loop over timestep (calc forces and propagate)
outputfile = open("1D_allosterysim.dat","w")
var_string = "#"

for i in range(N):
        var_string += "force[%d]   position[%d]   velocity[%d]   " %(i,i,i)
outputfile.write(var_string+"\n")

for it in range(Nstep):
	f = np.zeros(N)
	var_string = ""
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
		v[i] = c1*v[i] + np.random.normal()*(c2/m[i]) #O
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
	    
	outputfile.write(var_string+"\n")
	
# print output (coords, E,..)
outputfile.close()

# plot

data = np.loadtxt("1D_allosterysim.dat")

plt.hist(data[:,1],histtype='barstacked',label='Particle 1',bins=500)
plt.hist(data[:,4],histtype='step',label='Particle 2',bins=500)
plt.hist(data[:,7],histtype='barstacked',label='Particle 3',bins=500)
plt.hist(data[:,10],histtype='barstacked',label='Particle 4',bins=500)
plt.hist(data[:,13],histtype='step',label='Particle 5',bins=500)
plt.hist(data[:,16],histtype='barstacked',label='Particle 6',bins=500)
plt.hist(data[:,19],histtype='barstacked',label='Particle 7',bins=500)
plt.hist(data[:,22],histtype='step',label='Particle 8',bins=500)
plt.hist(data[:,25],histtype='barstacked',label='Particle 9',bins=500)

plt.xlabel("Positions",size=14)
plt.ylabel("Counts",size=14)
plt.legend(loc='upper left')

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

plt.savefig("9N_forceshist.pdf",dpi=600,transparant=True)
plt.close()

plt.hist(data[500:,2], histtype='step',label='Particle 1',bins=500)
plt.hist(data[500:,5], histtype='step',label='Particle 2',bins=500)
plt.hist(data[500:,8], histtype='step',label='Particle 3',bins=500)
plt.hist(data[500:,11], histtype='step',label='Particle 4',bins=500)
plt.hist(data[500:,14], histtype='step',label='Particle 5',bins=500)
plt.hist(data[500:,17], histtype='step',label='Particle 6',bins=500)
plt.hist(data[500:,20], histtype='step',label='Particle 7',bins=500)
plt.hist(data[500:,23], histtype='step',label='Particle 8',bins=500)
plt.hist(data[500:,26], histtype='step',label='Particle 9',bins=500)

plt.xlabel("Velocities",size=14)
plt.ylabel("Counts",size=14)
plt.legend(loc='upper left')

plt.savefig("9N_velocitieshist.pdf",dpi=600,transparant=True)
plt.close()


t = np.array(range(Nstep))*2

plt.plot(t,data[:,1],label='Particle 1')
plt.plot(t,data[:,4],label='Particle 2')
plt.plot(t,data[:,7],label='Particle 3')
plt.plot(t,data[:,10],label='Particle 4')
plt.plot(t,data[:,13],label='Particle 5')
plt.plot(t,data[:,16],label='Particle 6')
plt.plot(t,data[:,19],label='Particle 7')
plt.plot(t,data[:,22],label='Particle 8')
plt.plot(t,data[:,25],label='Particle 9')

plt.xlabel("Time (fs)",size=14)
plt.ylabel("Positions",size=14)
plt.legend(loc='upper left')
plt.savefig("9N_positionvstime.pdf",dpi=600,transparant=True)
plt.close()

t = np.array(range(Nstep-500))*2

plt.plot(t,data[500:,2],label='Particle 1')
plt.plot(t,data[500:,5],label='Particle 2')
plt.plot(t,data[500:,8],label='Particle 3')
plt.plot(t,data[500:,11],label='Particle 4')
plt.plot(t,data[500:,14],label='Particle 5')
plt.plot(t,data[500:,17],label='Particle 6')
plt.plot(t,data[500:,20],label='Particle 7')
plt.plot(t,data[500:,23],label='Particle 8')
plt.plot(t,data[500:,26],label='Particle 9')

plt.xlabel("Time (fs)",size=14)
plt.ylabel("Velocity",size=14)
plt.legend(loc='upper left')
plt.savefig("9N_velocitiesvstime.pdf",dpi=600,transparant=True)
plt.close()
