# Error Testing: Solution Accuracy vs. Analytical Soln
import solve as s
import plotting as pl
import potentials_1d as p
import pandas as pd
import numpy as np
import analytical as alt
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt

#----------------#
# Get User Input #
#----------------#
print(' Inft - 1\n SHO - 2\n')
type = int(input("Pick well type: "))


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

E_A = []
pct = []
if type ==1:
    U, x, v = p.infWell(N)

    #Numerical Computation:
    Ke = Cse*s.Laplacian(x)
    H = Ke + U
    E_N, psi = s.solution(H, N, x[0], x[N-1])

    #Analytical Computation:
    for n in range(len(E_N)):
        E_A.append(alt.infAlt(n))

    for i in range(len(E_N)):
        error = abs((E_A[i] - E_N[i])/E_A[i])
        pct.append(error*100)


elif type ==2:
    U, x, v = p.shopot(N)

    ##Numerical Computation:
    Ke = Cse*s.Laplacian(x)
    H = Ke + U
    E_N, psi = s.solution(H, N, x[0], x[N-1])

    #Analytical Computation:
    for n in range(len(E_N)):
        E_A.append(alt.shoAlt(n))

    for i in range(len(E_N)):
        error = abs((E_A[i] - E_N[i])/E_A[i])
        pct.append(error*100)

print(pct)

#---------------#
# Plot Solution #
#---------------#
plt.scatter(range(len(E_N)), pct, c="b", marker='x',
            label="Analytical")
#plt.scatter(range(len(E_N)), E_N, c="g", marker='o',
            #label="Numerical")

plt.xlabel("Energy Levels n")
plt.ylabel("pct error")
plt.legend(loc=2)
plt.show()
