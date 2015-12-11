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
data5000 = []
i = 0
# Fetching data from first calculation
with open('oppgave_MC_probability5000.txt') as oppgave_d:
	next(oppgave_d)
	for line in oppgave_d:
		x.append(float(line.split()[0])/100)
		data5000.append(float(line.split()[1]))
		i += 1

data200 = []
# Fetching data from first calculation
with open('oppgave_MC_probability200.txt') as oppgave_d:
	next(oppgave_d)
	for line in oppgave_d:
		data200.append(float(line.split()[1]))

data200g = []
# Fetching data from first calculation
with open('oppgave_MC_probability200gauss.txt') as oppgave_d:
	next(oppgave_d)
	for line in oppgave_d:
		data200g.append(float(line.split()[1]))

data5000g = []
# Fetching data from first calculation
with open('oppgave_MC_probability5000gauss.txt') as oppgave_d:
	next(oppgave_d)
	for line in oppgave_d:
		data5000g.append(float(line.split()[1]))

plt.figure(1)
plt.rc('font', family='serif')
plt.plot(x,data200,x,data5000,x,data200g,x,data5000g,xanal,analytical001,xanal,analytical025)
plt.legend([ur'$t=0.01$',ur'$t=0.25$',ur'$t_{\mathrm{gauss}}=0.01$',ur'$t_{\mathrm{gauss}}=0.25$',ur'$t_{\mathrm{an.}}=0.25$',ur'$t_{\mathrm{an.}}=0.01$'])
plt.title(ur"MC-simulering av diffusjon, n=99",fontsize=15)
plt.xlabel(ur'Position $x$',fontsize=16)
plt.ylabel(ur'Diffusion $u(x,t)$',fontsize=16)
plt.savefig('MC_walk.eps')
plt.show()