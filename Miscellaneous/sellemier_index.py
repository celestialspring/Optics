# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 11:31:39 2026

@author: SM
"""

import numpy as np
import matplotlib.pyplot as plt 


B = np.array([0.6962, 0.4074, 0.89748])
C = np.array([4.6791E-3, 1.3512E-2, 97.934]) #in um**2

ind = 3
def refractive_index(wavelength, B,C, ind):
    n1 = 0
    for i in range(ind):
        n1 += (B[i]*wavelength**2/(wavelength**2-C[i])) 
    
    n = np.sqrt(1+n1)
    return n
lambdas = np.array(np.linspace(0.4,0.7,5))

n = refractive_index(lambdas, B, C, ind)

print(n)

plt.plot(lambdas, n, 'g-',linewidth=2, markersize=12)
plt.xlabel('wavelength (um)')
plt.ylabel('Refractive index')