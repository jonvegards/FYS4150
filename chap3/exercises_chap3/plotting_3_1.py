#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Python script for plotting exercise 3.1

from matplotlib.pylab import *

with open('out.txt') as out:
        datapkt = out.readlines()

h = []
punkt = []

for line in datapkt:
    line.split()
    h.append = line[0]
    punkt.append = line[0]

plot(h,punkt)
show()
