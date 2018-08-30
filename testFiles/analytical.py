# Anylitical Solutions
import numpy as np
import math

#Analytical Computation:
def infAlt(n):
    n = n+1
    eV = 1.602*10**(-19)                                #kj*m^2/s^2
    kgnms = eV*(1*10**(18))                            #kj*nm^2/s^2
    hbar = (6.582119514e-16)*(kgnms)                  # Reduced Plancks  #kj*nm^2/s
    hbar2 = hbar**2
    me = 9.109e-31                          # Electon mass kg
    L2 = 0.1**2
    E = ((n**2*hbar2*np.pi**2)/(2*me*L2))
    En = E/(kgnms)
    return(En)

def shoAlt(n):
    hbar = 6.582119514e-16
    eV = 1.602*10**(-19)
    kgnms = eV*(1*10**(18))
    Vo = 400
    me = 9.109e-31
    L2 = 0.1**2
    w2 = (2*Vo)/(me*(L2))*(kgnms)
    w = math.sqrt(w2)
    E = hbar*w*(n + 0.5) - 400
    return(E)
