# -*- coding: utf-8 -*-
"""


@author: SM
"""
import numpy as np 
import tools
import matplotlib.pyplot as plt

class OpticalElements():
    
    def __init__(self, **kwargs):
        self.params = kwargs
        
    
    def rectslit(self, x_grid,y_grid, slitwidthx, slitwidthy):
        x_var = x_grid/(slitwidthx)
        y_var = y_grid/(slitwidthy)
        x_rect = OpticalElements.rectfunc(x_var)
        y_rect = OpticalElements.rectfunc(y_var)
        
        return x_rect, y_rect 
    
    @staticmethod
    def rectfunc(x):
        if isinstance(x, np.ndarray):  
          output = np.zeros_like(x, dtype=float)
          output[np.abs(x) == 0.5] = 0.5
          output[np.abs(x) < 0.5] = 1.0
          output[np.abs(x) > 0.5] = 0

        return output
    
    @staticmethod
    def circfunc(r):
        if isinstance(r, np.ndarray):
            output = np.zeros_like(r, dtype=float)
            output[r == 1] = 0.5
            output[r < 1 ] = 1
        
        return output
    
    def circslit(self, x_grid, y_grid, slit_diameter):
        r = np.sqrt(x_grid**2+y_grid**2)
        circ_var = (2*r)/slit_diameter
        circ_slit = OpticalElements.circfunc(circ_var)
        
        return circ_slit
    
if __name__ == '__main__':
    x = np.linspace(0,8E-3,100) -4E-3
    y = np.linspace(0,8E-3,100) -4E-3

    X,Y = np.meshgrid(x,y)
    clss = Opticalelements()
    xr, yr = clss.rectslit(X,Y, 100E-5)
    mask = xr*yr
    plt.figure()
    plt.imshow(mask, extent=[np.min(x), np.max(x), np.min(y), np.max(y)])
    plt.colorbar()
    plt.show()