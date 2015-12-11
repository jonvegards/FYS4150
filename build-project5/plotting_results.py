# -*- coding: utf8 -*-

import matplotlib.pyplot as plt
from numpy import array, linspace, pi, sin, exp

# Analytical solution
def AnalyticalSolution(x,t,N):
	u = linspace(0,1,len(x))
	tau = pi*x
	for i in range(0,len(x)):
		temp = 0
		for n in range(1,int(N)):
			temp2 = n*pi
			temp += -2*sin(n*pi*x[i])*exp(-(temp2)**2*t)/temp2
		u[i] = temp + 1 - x[i]
	return u

M=100
xanal = linspace(0,1,M)
tanal = linspace(0,1,M)
N = 1000

analytical001 = AnalyticalSolution(xanal,0.01, N)
analytical025 = AnalyticalSolution(xanal,0.25, N)

x = []
V1_implicit = []
V2_implicit = []

# Fetching data from first calculation
with open('oppgave_d_explicit_alfa.txt') as oppgave_d:
	next(oppgave_d)
	for line in oppgave_d:
		x.append(float(line.split()[0]))
		V1_implicit.append(float(line.split()[1])) 	# Short time
		V2_implicit.append(float(line.split()[2]))	# Steady state (almost)

x_exp = []
explicit1 = []
explicit2 = []

# Fetching data from first calculation
with open('oppgave_d_implicit_alfa.txt') as oppgave_b:
	next(oppgave_b)
	for line in oppgave_b:
		x_exp.append(float(line.split()[0]))
		explicit1.append(float(line.split()[1])) 	# Short time
		explicit2.append(float(line.split()[2]))	# Steady state (almost)

xCN = []
CN1 = []
CN2 = []

# Fetching data from first calculation
with open('oppgave_d_CN_alfa.txt') as oppgave_b:
	next(oppgave_b)
	for line in oppgave_b:
		xCN.append(float(line.split()[0]))
		CN1.append(float(line.split()[1])) 	# Short time
		CN2.append(float(line.split()[2]))	# Steady state (almost)

plt.figure(1)
plt.rc('font', family='serif')
plt.plot(x,V1_implicit,x,V2_implicit,x_exp,explicit1,x_exp,explicit2,xCN,CN1,xCN,CN2,xanal, analytical001,xanal, analytical025)
plt.legend(['Eksplisitt','Eksplisitt','Implisitt','Implisitt','Crank-Nicolson','Crank-Nicolson','t=0.01','t=0.25'])
plt.title(ur"Løsninger av diffusjonslikninga, n=99, $\Delta x = 1/10$, $\Delta t = \Delta x^2/2$",fontsize=15)
plt.xlabel(ur'Position $x$',fontsize=16)
plt.ylabel(ur'Diffusion $u(x,t)$',fontsize=16)
plt.savefig('oppgave_d_alfa.eps')

#plt.figure(2)
#plt.plot(xanal, analytical0,xanal, analytical05,xanal, analytical1)
#plt.legend(['t=0.01','t=0.1','t=1'])
plt.show()


