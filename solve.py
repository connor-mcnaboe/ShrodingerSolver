#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  File: solve.py
#  Created: 04/24/18
#  Copyright 2018 Connor <mithrandir@hodor>
#
#------------------------------------------------------------------------------#
# Solves  Shrodinger Equation using numpy linalg package                       #
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

import numpy as np
import math
import numpy.linalg as la

#-----------------#
# Simpsons Method #
#-----------------#

def simpsons(psi, a, b):
    #1-D Integration scheme: Simpsons rule.
    nS = len(psi)
    if nS % 2 == 1:
        sC = 2*np.ones(nS)
        sC[2:2:nS-1] = 4
        sC[1] = 1
        sC[nS-1] = 1
        h = (b-a)/(nS-1)
        sum = np.sum((h/3) * psi * sC)
    else:
        sum = "Not odd"
    return(sum)

#-----------#
# Laplacian #
#-----------#

def Laplacian(x):
    h = x[1]-x[0] #Delta X, step size
    n = len(x)                 #Delta X squared
    M = -2*np.identity(n,'d')
    for i in range(1,n):
        M[i,i-1] = M[i-1,i] = 1
    return M/h**2               #Differntial Matrix in one dimension

#------------------------#
#    Solve Shrodinger    #
#------------------------#

def solution(H, N, xMin, xMax):
    eigVal, eigVec = la.eigh(H) #Linalg Hermitian Matrix Calculation
    flag = 0
    n = 1
    E = []                      # Eigenvalues (Energies) within the well
    while flag == 0:
        e = eigVal[n-1]
        if e > 0:
            flag = 1
        else:
            E.append(e)
        n += 1                  #Number of eigenstates inside well
        if n > len(eigVal):
            flag = 1

    if not E:
        flag = 0
        n = 1
        for i in range(5):
            e = eigVal[i]
            E.append(e)
            n += 1
    # Find/Normalize the Eigenstates
    print(E)
    psi = np.zeros((N, n))
    for i in range(n-1):
        eVn = eigVec[:,i]
        eVn[0] = 0                                      #unnormalized B.Cs
        eVn[N-1] = 0
        psi[:,i] = eVn
        a = simpsons(eVn*eVn, xMin, xMax)               #Normalization
        psi[:,i] = eVn/math.sqrt(a)
        check = simpsons(psi[:,i]*psi[:,i], xMin, xMax) # Check if psi normal
        if check != 1:
            print('Normalization check not 1.0: ' + str(check))
    return(E, psi)

#-------------------------#
# Probablity Distribution #
#-------------------------#

def Probablity(E, psi):
    n = len(E)
    N = len(psi)
    prob = np.zeros((N, n))
    for i in range(n):
        prob[:, i] = psi[:,i]*psi[:,i]
    return(prob)
