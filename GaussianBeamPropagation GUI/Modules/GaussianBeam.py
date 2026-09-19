# -*- coding: utf-8 -*-
"""
Gaussian Beam Propagation

@author: SM
"""

import numpy as np
import matplotlib.pyplot as plt 

class GaussianBeamPropagation():
    
    def __init__(self, wavelength, waist=0, div_angle=0, rayleigh_range =0):
        
        self.wavelength = wavelength
        self.params = {'waist': waist, 'div_angle': div_angle, 'rayleigh_range' : rayleigh_range}
        non_zero_param = [name for name, val in self.params.items() if val !=0 ]
        print(non_zero_param)
        if non_zero_param[0] == 'waist':
            self.waist = self.params['waist']
            print(self.waist)
        elif non_zero_param[0] == 'div_angle':
            self.div_angle = self.params['div_angle']
        else:
            self.z_r = self.params['rayleigh_range']
    
        self.gaussian_param()
        
    def gaussian_param(self, n=1):
        
        lambda_0 = self.wavelength/n
        if self.waist:
            self.div_angle = lambda_0/(np.pi*self.waist) 
            self.z_r = self.waist/self.div_angle
        elif self.div_angle:
            self.z_r = lambda_0/(np.pi*self.div_angle**2)
            self.waist = self.z_r*self.div_angle
        else:
            self.waist = np.sqrt(lambda_0*self.z_r/(np.pi)) 
            self.div_angle = self.waist/self.z_r
        
        self.params['waist'] = self.waist
        self.params['div_angle'] = self.div_angle
        self.params['rayleigh_range'] = self.z_r
        
        
    def gaussian_field(self, x, z):
        n=1
        a_0 = 1
        lambda_0 = self.wavelength/n
        waist_z = self.waist* np.sqrt(1+(z/self.z_r)**2)
        radius_z = z * (1+(self.z_r/z)**2)
        phase = np.arctan(z/self.z_r)
        
        r = x
        k = 2*np.pi/(lambda_0)
        E_field = a_0 * (self.waist/waist_z) * np.exp(-r**2/(waist_z)**2 - 1j*k*z -1j*k*r**2/(2*radius_z) +1j*phase)
        
        Intensity = abs(E_field)**2
        
        return Intensity
    
    def simulation_plane(self):
        
        lx = 10 * self.waist
        lz = 3* self.z_r
        nz = 500 
        nx = 200
        x =  np.linspace(0, lx, nx)-lx/2
        z = np.linspace(0, lz, nz) 
        
        
        Z, X = np.meshgrid(z,x)
        
        return X,Z
    
    def plot(self, Intensity):
        X, Z = self.simulation_plane()
        norm_z = Z/self.z_r
        norm_x = X/self.waist
        plt.figure()
        plt.imshow(Intensity, cmap ='inferno', extent=[np.min(norm_z), np.max(norm_z), np.min(norm_x), np.max(norm_x)], aspect='auto')
        plt.colorbar(label='Intensity')
        plt.xlabel('$z/z_R$')
        plt.ylabel('$x/w_0$')
        plt.title("Gaussian beam propagation")

        info_text = f"$z_r$ = {round(self.z_r,2)} m, $w_0$ = {round(self.waist,4)} m"

        plt.text(0.25, 0.95, info_text, transform=plt.gca().transAxes, 
                 fontsize=10, verticalalignment='top', color='white')

if __name__ == '__main__':
    beam = GaussianBeamPropagation(632E-9, 0.5E-3)
    param = beam.params
    print(param)
    x, z = beam.simulation_plane()
    intensity = beam.gaussian_field(x,z)
    print(intensity.shape)
    beam.plot(intensity)
    
    