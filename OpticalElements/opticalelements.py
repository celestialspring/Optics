# -*- coding: utf-8 -*-
"""


@author: SM
"""
import numpy as np 
import tools
import matplotlib.pyplot as plt

class Opticalelements():
    
    def __init__(self, **kwargs):
        self.params = kwargs
        
    
    def rectslit(self, x_grid,y_grid, slitwidthx, slitwidthy):
        x_var = x_grid/(slitwidthx)
        y_var = y_grid/(slitwidthy)
        x_rect = Opticalelements.rectfunc(x_var)
        y_rect = Opticalelements.rectfunc(y_var)
        
        return x_rect, y_rect 
    
    @staticmethod
    def rectfunc(x):
        if isinstance(x, np.ndarray):  
          output = np.zeros_like(x, dtype=float)
          output[np.abs(x) == 0.5] = 0.5
          output[np.abs(x) < 0.5] = 1.0
          output[np.abs(x) > 0.5] = 0

        return output

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