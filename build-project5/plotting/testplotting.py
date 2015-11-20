import matplotlib.pyplot as plt
from numpy import array

x = []
V1 = []
V2 = []

# Fetching data from first calculation
with open('oppgave_d.txt') as oppgave_b:
	next(oppgave_b)
	for line in oppgave_b:
		x.append(float(line.split()[0]))
		V1.append(float(line.split()[1])) 	# Short time
		V2.append(float(line.split()[2]))	# Steady state (almost)


plt.plot(x,V1,x,V2)
plt.xlabel('Position')
plt.ylabel('Diffusion')
plt.show()