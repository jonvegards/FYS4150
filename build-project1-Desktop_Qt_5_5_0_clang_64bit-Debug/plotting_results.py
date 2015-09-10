#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
from numpy import log

# Script for plotting results from main.cpp

# Declaring lists
x = []
v = []
u = []
err = []

# Fetching data from first calculation
with open('oppgave_b_n_100000.txt') as oppgave_b:
	next(oppgave_b)
	for line in oppgave_b:
		x.append(float(line.split()[0]))
		v.append(float(line.split()[1])) 	# numerical
		u.append(float(line.split()[2]))	# exact
		err.append((float(line.split()[3])))
		#datapkt = oppgave_b.read().splitlines() #readlines()

# Fetching data from Armadillo calculation
x_arm = []
v_arm = []

with open('arm_solve_100000.txt') as oppgave_b:
        next(oppgave_b)
        for line in oppgave_b:
                x_arm.append(float(line.split()[0]))
                v_arm.append(float(line.split()[1])) 	# numerical
"""
plt.figure(1)
log, = plt.plot(x[2:-2], err[2:-2])
#plt.yscale('log')
plt.legend([log], ['Error'])
plt.title(r'Error plotted along $x$ for $n=100$')
plt.xlabel('$x_i$',size=17)
plt.ylabel('$\log_{10}(\epsilon)$',size=17)
plt.savefig("error_n_100.pdf")

"""
plt.figure(2)
ex,  = plt.plot(x,u, 'r', label='exact')
num, = plt.plot(x,v, 'b', label='numerical')
num_arm, = plt.plot(x_arm, v_arm, 'g', label='armadillo')
plt.legend([ex, num, num_arm], ['Exact', 'Numerical', 'Numerical w/ Armadillo'])
plt.xlabel('$x$',size=17)
plt.ylabel('$u(x)$',size=17)
plt.title(r'Results of calculation of $Ax = f$ for $n=100000$')
plt.savefig("d_n_100000.pdf")


"""
plt.figure(1)
num_arm, = plt.plot(x_arm, v_arm)
plt.legend([num_arm], ['Numerical'])
plt.title(r'Numerical solution for $n=10$')
"""
#plt.show()
