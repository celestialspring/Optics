# -*- coding: utf-8 -*-
"""
This script contains a class for calculating and plotting beam propagation for any given input field
using FFT-direct integration method.
It requires a config file for loading associated parameters, input field and grid sizes on which 
field is simulated.  

The script also contains a function to evaluate overlap between two modes.

@author: SM
"""

import tools
import mode_class
import numpy as np
import matplotlib.pyplot as plt 
import copy 
from scipy import integrate
import warnings

class PropagationMethods():
    
    def __init__(self, config_file:str, field:np.ndarray, X:np.ndarray, Y:np.ndarray):
        self.config_file = config_file
        self.field = field
        self.X = X
        self.Y = Y
        self.load_parameters()
    
    def load_parameters(self):
        
        self.param_dict = tools.read_configtxt(self.config_file)
        self.ef_amplitude = self.param_dict['field_amplitude']
        self.beam_waistr = self.param_dict['waist_radius']
        self.nm_order = self.param_dict['n_m']
        self.wavelength = self.param_dict['wavelength']
        self.points = self.param_dict['sample_points']
        self.grid_size = self.param_dict['grid_size']
        self.z_propagation = self.param_dict['propagation_distance']
        
    
    def fourier_propagator(self, kz,z):
        return np.exp(1j*kz*z)
    
    def real_space_propagator(self, k,R,z):
        self.h = (1/(2*np.pi)) * ((np.exp(1j*k*R ) * z) / R**2) * (1/R - (1j*k))
        return self.h
        
    def rayleigh_sommerfeld_DI(self):
        """
         
        Calculates diffraction integral using FFT
        Field is zero padded to 2N-1
        Grid points are mapped to 2N-1
        Ref: Applied optics Vol. 45, No. 6  20 February 2006   
        
        Returns
        -------
        np.ndarray 
            DESCRIPTION.

        """
        
        z = self.z_propagation
        nx = ny = int(self.points)
        lx, ly = self.grid_size
        #unit sampling length or pixel size
        dx = lx/nx
        dy = ly/ny
        
        #object plane N points
        XO = np.linspace(0,lx,nx) - lx/2
        YO = np.linspace(0,ly,ny) - ly/2
        
        #image plane N points
        xi = copy.deepcopy(XO)
        yi = copy.deepcopy(YO)
        
        # 2*n-1 points to account for padding real space propagator
        x_rh  = np.linspace(0,lx,2*nx-1) - lx/2
        y_rh  = np.linspace(0,ly,2*ny-1) - ly/2
        dX= np.zeros_like(x_rh) 
        dY= np.zeros_like(y_rh)
        
        for i in range(nx-1): #till n-2
            dX[i] = xi[0] - XO[nx-1-i]  
            dY[i] = yi[0] - YO[nx-1-i]  
        
        for j in range(nx-1,2*nx-1):
            dX[j] = xi[j-nx+1] - XO[0] 
            dY[j] = yi[j-nx+1] - YO[0] 
           
        DX,DY = np.meshgrid(dX,dY)

        R = np.sqrt(DX**2+DY**2+z**2)
        k = (2*np.pi)/self.wavelength
                    
        #padding input field
        U_pad = np.zeros((2*ny-1, 2*nx-1), dtype=complex)
        U_pad[0:ny, 0:nx] = self.field
        
        if np.isclose(z,0):
            self.S_shift = U_pad[0:ny-1, 0:nx-1]
            warnings.warn('Returning initial field, propagation length is zero !',UserWarning)
            return self.S_shift 
        
        H = self.real_space_propagator(k, R, z)
        S_fourier = np.fft.fftshift(np.fft.fft2(U_pad)) * np.fft.fftshift(np.fft.fft2(H))
        S = np.fft.ifft2(np.fft.ifftshift(S_fourier)) *dx * dy
        self.S_shift = S[ny-1:2*ny-1,nx-1:2*nx-1] # dx,dy=0 is at lower right somewhere near nx
        return self.S_shift
    
    def plot_field(self, plot_type: str):
        
        xsize = format(self.grid_size[0],'f')
        ysize = format(self.grid_size[1],'f')
        lx, ly = self.grid_size
        x_str = str(xsize)
        y_str = str(ysize)
        
        x_decimal = len(x_str.split('.')[1])
        y_decimal = len(y_str.split('.')[1])
       
        X = self.X*(10**x_decimal)
        Y = self.Y*(10**y_decimal)
        if plot_type == 'Amplitude':
          
            plt.figure()
            plt.imshow(abs(self.S_shift),extent=[np.min(X), np.max(X), np.min(Y), np.max(Y)])
            plt.colorbar(label='Amplitude')
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title(f'HG mode {int(self.nm_order[0]),int(self.nm_order[1])} field amplitude, z={self.z_propagation} m')
            plt.show()
            
        if plot_type == 'Intensity':
            plt.figure()
            field_intensity = (abs(self.S_shift))**2
            plt.imshow(field_intensity,extent=[np.min(X), np.max(X), np.min(Y), np.max(Y)])
            plt.colorbar(label='Intensity')
            plt.xlabel('x')
            plt.ylabel('y')
            plt.title(f'HG mode {int(self.nm_order[0]),int(self.nm_order[1])} field intensity, z={self.z_propagation} m')
            plt.show()
        
        if plot_type == 'Intensity xprofile':
            
            field_intensity = (abs(self.S_shift))**2
            nxi, nyi = field_intensity.shape
            mid_y = nyi // 2
            inter_x =0.5 * (X[0, :-1] + X[0, 1:])#grid points at vertices, Intensity at mid_points
            plt.figure()
            plt.plot(X[0,:], field_intensity[mid_y, :])
            plt.xlabel('x')
            plt.ylabel('Intensity')
            plt.title(f'HG mode {int(self.nm_order[0]),int(self.nm_order[1])} intensity profile, z={self.z_propagation} m')
            plt.show()

def two_mode_overlap(mode1, mode2, xgrid, ygrid):
    """
    Calculates overlap of two modes 

    Returns
    -------
    overlap : float
       scalar value of the overlap (0-1)

    """
    
    if mode1.shape != mode2.shape:
        raise Exception('Sampling points must be same for the two modes !')
        
    # xgrid_inter =0.5 * (xgrid[0, :-1] + xgrid[0, 1:])
    # ygrid_inter =0.5 * (ygrid[:-1, 0] + ygrid[1:, 0])
    xgrid_inter = xgrid[0,:]
    ygrid_inter = ygrid[:,0]
    inner_int_y = integrate.trapezoid(mode1*np.conj(mode2),ygrid_inter, axis=0)
    inner_product_square = np.abs((integrate.trapezoid(inner_int_y,xgrid_inter, axis=0)))**2

    mode1_int_y = integrate.trapezoid(np.abs(mode1)**2,ygrid_inter, axis=0)
    mode1_power = (integrate.trapezoid(mode1_int_y,xgrid_inter, axis=0))
    
    mode2_int_y = integrate.trapezoid(np.abs(mode2)**2,ygrid_inter, axis=0)
    mode2_power = (integrate.trapezoid(mode2_int_y,xgrid_inter, axis=0))

    overlap = round(inner_product_square/(mode1_power*mode2_power),3)

    return overlap
    
if __name__ == '__main__':
    
    new_mode = mode_class.Modes('mode1_config.txt')
    X,Y = new_mode.generate_object_plane()
    field = new_mode.field_input
    mode_propagated = PropagationMethods('mode1_config.txt', field, X, Y)
    mode_propagated.rayleigh_sommerfeld_DI()
    mode_propagated.plot_field('Intensity xprofile')
    mode_propagated.plot_field('Intensity')