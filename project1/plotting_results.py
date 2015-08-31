#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pylab as plt

# Script for plotting results from main.cpp

# Declaring lists
x = []
v = []
u = []
err = []

# Opening file
with open('oppgave_b_n_1000.txt') as oppgave_b:
	next(oppgave_b)
	for line in oppgave_b:
		x.append(float(line.split()[0]))
		v.append(float(line.split()[1]))
		u.append(float(line.split()[2]))
		err.append(abs(float(line.split()[3])))
		#datapkt = oppgave_b.read().splitlines() #readlines()


plt.figure(1)
plt.plot(x, err)
plt.yscale('log')
#plt.xscale('log')

plt.figure(2)
plt.plot(x,u, 'r')

plt.plot(x,v, 'b')

plt.show()
