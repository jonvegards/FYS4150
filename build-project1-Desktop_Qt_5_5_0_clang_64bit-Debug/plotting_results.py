#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pylab as plt
from numpy import log

# Script for plotting results from main.cpp
"""
# Declaring lists
x = []
v = []
u = []
err = []

# Opening file
with open('oppgave_b_n_10000.txt') as oppgave_b:
	next(oppgave_b)
	for line in oppgave_b:
		x.append(float(line.split()[0]))
		v.append(float(line.split()[1])) 	# numerical
		u.append(float(line.split()[2]))	# exact
		err.append((float(line.split()[3])))
		#datapkt = oppgave_b.read().splitlines() #readlines()


plt.figure(1)
log, = plt.plot(x[2:-2], err[2:-2])
#plt.yscale('log')
plt.legend([log], ['Error'])
plt.title(r'Error plotted along $x$ for $n=100$')

plt.figure(2)
ex,  = plt.plot(x,u, 'r', label='exact')
num, = plt.plot(x,v, 'b', label='numerical')
plt.legend([ex, num], ['Exact', 'Numerical'])
plt.title(r'Results of calculation of $Ax = f$ for $n=100$') """

#############################################################################
# Plotting results from Armadillo run

x_arm = []
v_arm = []

with open('arm_solve_10.txt') as oppgave_b:
        next(oppgave_b)
        for line in oppgave_b:
                x_arm.append(float(line.split()[0]))
                v_arm.append(float(line.split()[1])) 	# numerical


plt.figure(1)
num_arm, = plt.plot(x_arm, v_arm)
plt.legend([num_arm], ['Numerical'])
plt.title(r'Numerical solution for $n=10$')

plt.show()
