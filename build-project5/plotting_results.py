import matplotlib.pyplot as plt
from numpy import array, linspace, pi, sin, exp

# # Analytical solution
# def AnalyticalSolution(x,t,N):
# 	temp = 0
# 	u = linspace(0,1,M)
# 	Lambda = 1
# 	d = 1
# 	tau = 2*pi*x
# 	for n in range(0,len(x)):
# 		for i in range(0,int(N)):
# 			temp += sin(i*x[n])
# 		u[n] = (temp/N + 1 - x[n])*exp(-Lambda*(t)**2)
# 		temp = 0
	
# 	return u

# M=1000
# xanal = linspace(0,1,M)
# tanal = linspace(0,1,M)
# N = 2000

# #analytical0 = AnalyticalSolution(xanal,0, N)
# analytical1 = AnalyticalSolution(xanal,1, N)
# print analytical1

x = []
V1_implicit = []
V2_implicit = []

# Fetching data from first calculation
with open('oppgave_d_explicit.txt') as oppgave_d:
	next(oppgave_d)
	for line in oppgave_d:
		x.append(float(line.split()[0]))
		V1_implicit.append(float(line.split()[1])) 	# Short time
		V2_implicit.append(float(line.split()[2]))	# Steady state (almost)

x_exp = []
explicit1 = []
explicit2 = []

# Fetching data from first calculation
with open('oppgave_d_implicit.txt') as oppgave_b:
	next(oppgave_b)
	for line in oppgave_b:
		x_exp.append(float(line.split()[0]))
		explicit1.append(float(line.split()[1])) 	# Short time
		explicit2.append(float(line.split()[2]))	# Steady state (almost)

xCN = []
CN1 = []
CN2 = []

# Fetching data from first calculation
with open('oppgave_d_CN.txt') as oppgave_b:
	next(oppgave_b)
	for line in oppgave_b:
		xCN.append(float(line.split()[0]))
		CN1.append(float(line.split()[1])) 	# Short time
		CN2.append(float(line.split()[2]))	# Steady state (almost)

plt.rc('font', family='serif')
plt.plot(x,V1_implicit,x,V2_implicit,x_exp,explicit1,x_exp,explicit2,xCN,CN1,xCN,CN2)#, xanal, analytical1)
plt.legend(['Eksplisitt','Eksplisitt','Implisitt','Implisitt','Crank-Nicolson','Crank-Nicolson'])#,'Analytical'])
plt.title(ur"Solutions of diffusion equation, n=99, $\Delta x = 1/10$, $\Delta t = \Delta x^2/2$",fontsize=15)
plt.xlabel(ur'Position $x$',fontsize=16)
plt.ylabel(ur'Diffusion $u(x,t)$',fontsize=16)
plt.savefig('oppgave_d.eps')
plt.show()