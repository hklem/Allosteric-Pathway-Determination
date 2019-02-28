import numpy as np
import matplotlib.pyplot as plt


m = 12. # amu
T = 400 # K
k = 0.00198715 

data = np.loadtxt("1D_allosterysim.dat")
v = np.arange(-1,1,0.001)
mbd = np.sqrt(m/(2*np.pi*k*T))*np.exp((-m*v**2)/(2*k*T))


v_var = np.var(data[500:,2])

mbd2 = np.sqrt(1./(2.*np.pi*v_var))*np.exp((-v**2)/(2.*v_var))

Temp = (v_var*m)/k
print(Temp)

plt.hist(data[500:,2],density=True,bins=50)
plt.plot(v,mbd)
plt.plot(v,mbd2)
plt.savefig("maxwell_dist.pdf",dpi=600)
plt.close()

