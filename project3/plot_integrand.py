# Python-script for plotting function to be integrated in project 3

import numpy as np
import matplotlib.pylab as plt

Limit = .5
n = 100

x1 = np.linspace(-Limit,Limit,n)
y1 = np.linspace(-Limit,Limit,n)
z1 = np.linspace(-Limit,Limit,n)
x2 = np.linspace(-Limit,Limit,n)
y2 = np.linspace(-Limit,Limit,n)
z2 = np.linspace(-Limit,Limit,n)

def Integrand(x1,y1,z1,x2,y2,z2):
	r1vec = np.array((x1,y1,z1),dtype=float)
	r2vec = np.array((x2,y2,z2),dtype=float)
	r1 = np.sqrt(r1vec[0,:]**2 + r1vec[1,:]**2 + r1vec[2,:]**2)
	r2 = np.sqrt(r2vec[0,:]**2 + r2vec[1,:]**2 + r2vec[2,:]**2)
	return np.exp(-2*2*(r1+r2)), r1

funkis, r1 = Integrand(x1,y1,z1,x2,y2,z2)

plt.plot(r1,funkis)
plt.show()