#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  File: main.py
#  Created: 04/24/18
#  Copyright 2018 Connor <mithrandir@hodor>
#
#-------------#
# Main Script #
#-------------#
import solve as s
import plotting as pl
import potentials_1d as p
import pandas as pd
import numpy as np
#import testing as tst
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt

#----------------#
# Get User Input #
#----------------#
print(' Inft - 1\n Square - 2\n Stepped - 3\n Sloped -4\n Lattice - 5\n SHO - 6\n')
type = int(input("Pick well type: "))
print(' Psi - 1\n Probablity -2\n')
ptype = int(input('Pick a plot type: '))
level = int(input('Choose energy level: '))
#Universial Constants#

hbar =  1.055e-34                       # Reduced Plancks  J.s
Lse = 1e-9                              # m-> nm
Ese = 1.6e-19                           #J -> eV

#Particle-Specific Constants# (Idea is to make file (JS?) stroring the data)

#Electron
me = 9.109e-31                          # Electon mass kg
e = 1.602e-19                           # One columb
Cse = -(hbar**2/(2*me))/((Lse**2*Ese))  #Shrodinger Constant

N = 1001 #Number of disctretaztions

#-----------------#
# Build Matricies #
#-----------------#

#Sloped Well
if type ==1:
    U, x, v = p.infWell(N)
elif type ==2:
    U, x, v = p.squareWell(N)
elif type==3:
    U, x, v = p.stepWell(N)
elif type==4:
    U, x, v = p.slopedWell(N)
elif type==5:
    U, x, v = p.lattice(N)
elif type ==6:
    U, x, v = p.shopot(N)


#------Opertations that occur regardless of input----#

Ke = Cse*s.Laplacian(x)
H = Ke + U
E, psi = s.solution(H, N, x[0], x[N-1])
prob = s.Probablity(E, psi)

#---------------#
# Plot Solution #
#---------------#

edf = pd.DataFrame(E)
print(E[level - 1])

xMin = np.min(x)
xMax = np.max(x)
mtitle, ylabel = pl.plot_labels(type, ptype)
host, par1 = pl.plot_setup(mtitle, ylabel, xMin, xMax)

if ptype == 1:
    yd = psi[:,level-1]
else:
    yd = prob[:,level-1]

p1, = host.plot(x, yd, label=ylabel)
y1min = np.min(yd) - 5
y1max = np.max(yd) + 5
y2min = np.min(v) - 10
y2max = np.max(v) + 10
if type == 1:
    p2, = par1.plot(x, np.zeros(N))
else:
    p2, = par1.plot(x, v, label="Potential")

pl.finish_plot(host, par1, p1, p2, y1min, y1max, y2min, y2max)
