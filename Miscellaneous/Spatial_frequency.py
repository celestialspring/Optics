# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 16:21:06 2026

@author: SM
"""

import numpy as np
import matplotlib.pyplot as plt 
'Grid'
N=500
L = 20
xoverlambda = np.linspace(-10, 10, N)
zoverlambda = np.linspace(0,20,N)
dzlambda = dxlambda = L/N
ZLam, XLam = np.meshgrid( zoverlambda, xoverlambda)


'fourier coeff'
fxlambda_f = 0.5
fzlambda_f = np.lib.scimath.sqrt(1-fxlambda_f**2) #handles imaginary
dfxlambda = 1/dxlambda

fzlambda = fxlambda = np.fft.fftshift(np.fft.fftfreq(N, dxlambda))
'function'
E_0 = 1
E = E_0 *np.exp(1j*2*np.pi*(fxlambda_f*XLam+fzlambda_f*ZLam))
ft_E = np.fft.fftshift(np.fft.fft(E[:,0]))

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 8))

axes[0,1].imshow(np.real(E),extent=[np.min(zoverlambda),np.max(zoverlambda),np.min(xoverlambda),np.max(xoverlambda)])
axes[0,1].set_title('Real space plane wave')
axes[0,0].plot(np.real(E[:,0]), xoverlambda)
axes[0,0].set_box_aspect(1)
ticky= np.linspace(-10,10,11)
axes[0,0].set_yticks(ticky)
axes[0,0].set_title('x/$\lambda$ vs Re(E(x,z=0))')
axes[1,0].plot(zoverlambda,np.real(E[0,:]))
axes[1,0].set_box_aspect(1)
axes[1,0].set_title('z/$\lambda$ vs Re(E(x=0,z))')
axes[1,1].plot((ft_E)/np.max(abs(ft_E)),fxlambda)
axes[1,1].set_box_aspect(1.5)
axes[1,1].set_title('$f_{x}$$\lambda$ vs FT(E($f_{x}$,z=0))')
fig.suptitle('Spatial Frequency Explorer Overview', fontsize=14)