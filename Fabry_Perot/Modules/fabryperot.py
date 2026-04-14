# -*- coding: utf-8 -*-
"""
A class that contains relevant equations to calculate FP parameters

@author: SM
"""
import numpy as np 

class FabryPerot():
    
    def __init__(self, wavelength, r, cavity_length):
        self.wavelength = wavelength
        self.r = r
        self.d = cavity_length
        
    def coeff_Of_finesse(self):
        F = 4*self.r**2/(1-self.r**2)**2
        return F
    
    def finesse(self):
        """
        
        It is also the ratio between the separation of peaks to FWHM of the peaks.
        In other words, it shows how well the cavity is tuned or emphasizes resonance effect of the cavity 
        Returns
        -------
        Fi : TYPE
            DESCRIPTION.

        """
        
        Fi = np.pi*self.r/(1-self.r**2)
        return Fi
    
    def transmittance(self, F):
        
        delta = 2*2*np.pi*self.d/(self.wavelength)
        T = 1/(1+F*np.sin(delta/2)**2)
        return T
    
        