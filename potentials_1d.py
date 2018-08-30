#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  File: potentials.py
#  Created: 04/24/18
#  Copyright 2018 Connor <mithrandir@hodor>
#
#
#------------------------------------------------------------------------------#
# Initalizes Potential Energy Matricies                                        #
# Function Defined:                                                            #
#   simpsons(psi, a, b) - Numerical integration used to normalize Psi          #
#       Returns: sum (int) squared normalization coefficient                   #
#   Laplacian(x) - Creates diagonalized Second derivative matrix               #
#       Returns: M (ndarray) d^2/dx^2                                          #
#   solution() - Solves for eigenvalues and eigenvectors                       #
#       Returns: E, psi(narray,ndarray) normalized Eigenstates & well energies #
#   probability() - Creates Probablity distribution for each States            #
#       Returns: prob (ndarray) Matrix, each column represents diff state      #
#------------------------------------------------------------------------------#


#=================Notes==============================#
# Some list/dict input storing system specifications #
# would be the imput, each indice corresponding to a #
# adjustable variable. Hardcoded optins temp(Input)  #
#----------------------------------------------------#

import numpy as np
import math

#---------------------#
# Initalize Constants #
#---------------------#

#Universial Constants#

hbar =  1.055e-34                       # Reduced Plancks  J.s
Lse = 1e-9                              # m-> nm
Ese = 1.6e-19                           #J -> eV

#Particle-Specific Constants# (Idea is to make file (JS?) stroring the data)

#Electron
me = 9.109e-31                          # Electon mass kg
e = 1.602e-19                           # One columb
Cse = -(hbar**2/(2*me))/((Lse**2*Ese))  #Shrodinger Constant

#-------------------------#
# Infinite Well Potential #
#-------------------------#
def infWell(N):
    v0 = 0                                  # Square Well Depth: defualt = -400 eV
    v = np.zeros(N)                         # Potential array
    U = np.zeros([N, N])                    # Potnetial matix
    xMin = -0.05                            # Left boundry of well: defualt = -0.05nm
    xMax = 0.05                             # Right Boundry of well: Defualt = 0.05nm
    xW2 = 0.025                             # 1/2 width of the well: defualt = 0.025nm
    x = np.linspace(xMin, xMax, num=N)      # x values


    for i in range(N-1):                    #Square well array setup
        if abs(x[i]) <= xW2:
            v[i] = v0

    for i in range(N-1):                    #Square well matrix
        for j in range(N-1):
            if i == j:
                U[i,j] = v[j]
    return(U, x, v)

#-----------------------#
# Square Well Potential #
#-----------------------#
def squareWell(N):
    v0 = -400.00                            # Square Well Depth: defualt = -400 eV
    v = np.zeros(N)                         # Potential array
    U = np.zeros([N, N])                    # Potnetial matix
    xMin = -0.1                             # Left boundry of well: defualt = -0.1nm
    xMax = 0.1                              # Right Boundry of well: Defualt = 0.1nm
    xW2 = 0.05                              # 1/2 width of the well: defualt = 0.05nm
    x = np.linspace(xMin, xMax, num=N)      # x values


    for i in range(N-1):                    #Square well array setup
        if abs(x[i]) <= xW2:
            v[i] = v0

    for i in range(N-1):                    #Square well matrix
        for j in range(N-1):
            if i == j:
                U[i,j] = v[j]
    return(U, x, v)

#---------------------#
# Step Well Potential #
#---------------------#
def stepWell(N):
    xMin = -0.15            # default = -0.15 nm
    xMax = 0.15             # default = 0.15 nm
    x1 = 0.1                # Total width of well: default = 0.1 nm
    x2 = 0.06               # Width of LHS well: default = 0.06 nm
    U1 = -400               # Depth of LHS well: default = -400 eV
    U2 = -250               # Depth of RHS well (eV):default = -250 eV
    v = np.zeros(N)         # Potential array
    U = np.zeros([N, N])    # Potnetial matix

    x = np.linspace(xMin,xMax, num=N)

    for i in range(N-1):
        if x[i] >= -x1/2 and x[i] < -x1/2 + x2:
             v[i] = U1
        elif x[i] >= x1/2 - x2 and x[i] < x1/2:
            v[i] = U2
        else:
            v[i] = 0

    for i in range(N-1):                    #Diagonalzing
        for j in range(N-1):
            if i == j:
                U[i,j] = v[j]
    return(U, x, v)


#-----------------------#
# Sloped Well Potential #
#-----------------------#

def slopedWell(N):
    xMin = -0.1      # default value = -0.1 nm
    xMax = 0.1       # default value = + 0.1 nm
    U1 = -1200       # Depth of LHS well: default = -1200 eV;
    U2 = -200        # Depth of RHS well: default = -200 eV;
    x1 = 0.07        # 1/2 width of well: default = 0.05 nm;
    v = np.zeros(N)         # Potential array
    U = np.zeros([N, N])    # Potnetial matix

    #Set Linear eq
    x = np.linspace(xMin,xMax, num=N)
    intercept = (U1+U2)/2         #Middle of well
    slope = (U2-U1)/(2*x1)
    for i in range(N-1):
        if abs(x[i])<= x1:
            v[i] = slope * x[i] + intercept

    for i in range(N-1):                    #Diagonalzing
        for j in range(N-1):
            if i == j:
                U[i,j] = v[j]
    return(U, x, v)


def lattice(N):
    num = N
    wellNum = 12      # number of wells
    U1 = -350
    x1 = 0.05
    x2 = 0.075
    xEnd = 0.05
    wellDepth = U1*np.ones(wellNum)   # depth of wells
    wellWidth = x1*np.ones(wellNum)
    wellSeparation = x2*np.ones(wellNum-1)
    wellcenter = np.zeros(wellNum)
    xMin = 0
    xMax = 2*xEnd + np.sum(wellSeparation) + 0.5*(wellWidth[1]+wellWidth[wellNum-1])
    x = np.linspace(xMin,xMax,num)
    dx = (xMax-xMin)/(num-1)
    v = np.zeros(N)         # Potential array
    U = np.zeros([N, N])
    wellcenter[1] = xMin+xEnd+wellWidth[1]/2

    for i in range(2, wellNum-1):
        wellcenter[i] = wellcenter[i-1] + wellSeparation[i-1]


    for i in range(wellNum-1):
        for j in range(num-1):
            if abs(x[j]-wellcenter[i]) <= wellWidth[i]/2:
                v[j] = wellDepth[i]

    for i in range(N-1):                    #Diagonalzing
        for j in range(N-1):
            if i == j:
                U[i,j] = v[j]
    return(U, x, v)

def shopot(N):
    xMin = -0.2                  # default = -0.0 nm
    xMax = 0.2                   # default = +0.2 nm
    x1 = 0.2                     # width default = 0.2 nm
    U1 = -400                    # well depth default = -400 eV
    v = np.zeros(N)              # Potential array
    U = np.zeros([N, N])         # Potnetial matix

    x = np.linspace(xMin,xMax, num=N)
    for i in range(N):
        if abs(x[i])<= x1/2:
            v[i] = -(4*U1/(x1*x1))*x[i]**2+U1

    for i in range(N-1):                    #Diagonalzing
        for j in range(N-1):
            if i == j:
                U[i,j] = v[j]
    return(U, x, v)
