# -*- coding: utf-8 -*-
"""
This script contains a class used to generate and plot field of a beam mode.
It can generate HG mode upto 3 orders.
Parameters are taken from a config file given as input.
@author: SM
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        """
import tools 
import numpy as np
import matplotlib.pyplot as plt 
import warnings

class Modes():
    
    def __init__(self, config:str):
        
        self.config_file = config
        self.load_parameters()
        self.field_input = self.select_mode()
    
    def load_parameters(self):
        
        self.param_dict = tools.read_configtxt(self.config_file)
        self.ef_amplitude = self.param_dict['field_amplitude']
        self.beam_waistr = self.param_dict['waist_radius']
        self.nm_order = self.param_dict['n_m']
        self.wavelength = self.param_dict['wavelength']
        self.points = self.param_dict['sample_points']
        self.grid_size = self.param_dict['grid_size']
        self.z_propagation = self.param_dict['propagation_distance']
        self.rayleigh_length = (np.pi * self.beam_waistr**2)/self.wavelength
        self.radius = self.beam_waistr*np.sqrt(1+(self.z_propagation/self.rayleigh_length)**2)
        

    def generate_object_plane(self):
        """
        Generates a grid of points in x-y plane

        Returns
        -------
        X : X-coordinate of a grid
         
        Y : Y-coordinate of a grid
            
        """
        nx = ny = int(self.points)
        lx, ly = self.grid_size
        
        x = np.linspace(0,lx,nx) - lx/2
        y = np.linspace(0,ly,ny) - ly/2
        
        X,Y = np.meshgrid(x,y)
        return X,Y
    
    
    def grid_kcomponents(self):
        """
        Calculates possible k components corresponding to the grid size and 
        sampling parameters

        Returns
        -------
        k, kX,kY,kZ : np.ndarray
        """
                
        nx = ny = int(self.points)
        lx, ly = self.grid_size
        
        #unit sampling length or pixel size
        if lx==ly:
            dx = dy = lx/nx
            #n/dx*nx, gives sampling frequencies indexed by n and centre it with fftshift
            kbar = np.fft.fftshift(np.fft.fftfreq(nx, dx))
            kbar_x = kbar_y = kbar            
        else:
            dx = lx/nx
            dy = ly/ny
            kbar_x = np.fft.fftshift(np.fft.fftfreq(nx, dx))
            kbar_y = np.fft.fftshift(np.fft.fftfreq(ny, dy))
        
        #wave_vectors 
        kx = 2*np.pi*kbar_x
        ky = 2*np.pi*kbar_y
        kX, kY = np.meshgrid(kx, ky)
        k = (2*np.pi)/self.wavelength
        kZ = np.sqrt(k**2-kX**2-kY**2, dtype='complex')
        
        return k, kX, kY, kZ
    
    def sampling_check(self):
        """
        Check if sampling/grid is sufficiently large
        f_s > 2*f_max (Nyquist), spatial units [L-1]
        f_max < 1/2*dx
        resolution check dx< w0/8
        grid size check L> 4*w0
        Returns
        -------
        None.
        
        """
        nx = ny = int(self.points)
        lx, ly = self.grid_size
        k, kX, kY, kZ = self.grid_kcomponents()
        if not np.isclose(self.z_propagation, 0):
            theta = np.arctan(self.radius/self.z_propagation) #gaussian divergence
        else:
            theta = 0
            warnings.warn("propagation distance too low", UserWarning)
        kx = k*theta
        f_max = kx/(2*np.pi) 
        if lx == ly:
           dx = lx/nx
           fs = 1/dx
           if dx > self.beam_waistr/8:
               warnings.warn(f"Weak spatial resolution dx = {dx} ! increase sampling points", UserWarning)
           if f_max <= fs/2 and lx >= (4*self.beam_waistr):
               pass
           else:
                warnings.warn("Weak sampling !", UserWarning)
                if f_max > fs/2:
                    print(f'f_max > f_s/2, fmax = {f_max}, fs ={fs}')
                if lx < (4*self.beam_waistr):
                    print(f'grid size < 4*beam_waist, 4*w0 ={4*self.beam_waistr}')      
        else:
            dx = lx/nx
            dy = ly/ny
            fsx = 1/dx
            fsy = 1/dy
            
            if min(dx,dy) >  self.beam_waistr/8:
                warnings.warn("Weak spatial resolution ! increase sampling points", UserWarning)
            if min(f_max,f_max) <= min(fsx/2, fsy/2)  and min(lx,ly) >= 4*self.beam_waistr:
                pass
            else:
                warnings.warn("Weak sampling !", UserWarning)
                if f_max > fsx/2 or f_max > fsy/2:
                    print(f'f_max > f_s/2, fmax = {f_max}, fsx ={fsx}, fsy ={fsy}')
                if lx < (4*self.beam_waistr) or ly < (4*self.beam_waistr):
                    print(f'grid size < 4*beam_waist, 4*w0 ={4*self.beam_waistr}')
        
    def select_mode(self):
        """
        Selects a given mode type specified in config file
        
        Returns
        -------
        None.

        """
        if self.param_dict['mode_type'] == 'Hermite-Gaussian':
            self.X, self.Y = self.generate_object_plane()
            z = 0 
            self.sampling_check()
            return self.mode_hg(self.X,self.Y,z)
            
    
    def hermite_function(self, order:int, x):
        """

        Parameters
        ----------
        order : int
            order of the hermite polynomial.
        x : TYPE
            Hermite polynomial is a function of x; H(x).

        Raises
        ------
        Exception
            Cannot go beyond n=3.

        Returns
        -------
        h : float
            returns calculated polynomial.

        """
        n = order
        if n == 0:
            h = 1
        elif n == 1:
            h = 2*x
        elif n == 2:
            h = 4*x**2 - 2
        elif n == 3:
            h = 8*x**3 - 12*x
        else:
            raise Exception('Order cannot exceed 3')
        return h     
    
    def mode_hg(self,X,Y,z:float):
        """
        Field corresponding to HG mode is calculated

        Returns
        -------
        np.ndarray
            DESCRIPTION.

        """
        
        X0 = Y0 = 0
        n = self.nm_order[0]
        m = self.nm_order[1]
        
        'hermite polynomials'
        x_hermite = (np.sqrt(2) * X )/ self.beam_waistr
        y_hermite = (np.sqrt(2) * Y )/ self.beam_waistr
        hn = self.hermite_function(n, x_hermite)
        hm = self.hermite_function(m, y_hermite)
        
        r = np.sqrt((X-X0)**2+(Y-Y0)**2)
        R = z*(1+(z/self.rayleigh_length)**2) #wavefront curvature
        k = (2*np.pi)/self.wavelength
        
        "Electric field distribution"
        self.field_eq = self.ef_amplitude*(self.rayleigh_length/self.beam_waistr)*hn*hm\
        * np.exp(-r**2/self.beam_waistr**2)\
        * np.exp(1j*k*z-(1+n+m)*np.arctan(z/self.rayleigh_length)+ (k*r**2)/2*R)
        return self.field_eq
    
    def plot_field(self, plot_type: str):
        
        xsize = format(self.grid_size[0],'f')
        ysize = format(self.grid_size[1],'f')
        
        x_str = str(xsize)
        y_str = str(ysize)
        
        x_decimal = len(x_str.split('.')[1])
        y_decimal = len(y_str.split('.')[1])
        
        X = self.X*(10**x_decimal)
        Y = self.Y*(10**y_decimal)
        
        if plot_type == 'Amplitude':     
            plt.figure()
            plt.imshow(abs(self.field_eq),extent=[np.min(X), np.max(X), np.min(Y), np.max(Y)])
            plt.colorbar(label='Amplitude')
            plt.xlabel('x ')
            plt.ylabel('y ')
            plt.title('Hermite-Gaussian field amplitude')
            plt.show()
            
        if plot_type == 'Intensity':
            plt.figure()
            field_intensity = abs(self.field_eq)**2
            plt.imshow(field_intensity,extent=[np.min(X), np.max(X), np.min(Y), np.max(Y)])
            plt.colorbar(label='Intensity')
            plt.xlabel('x ')
            plt.ylabel('y ')
            plt.title('Hermite-Gaussian field intensity')
            plt.show()
            
        if plot_type == 'Intensity xprofile':
            
            field_intensity = (abs(self.field_eq))**2
            nxi, nyi = field_intensity.shape
            mid_y = nyi // 2
            plt.figure()
            plt.plot(X[0,:], field_intensity[mid_y, :])
            plt.xlabel('x')
            plt.ylabel('Intensity')
            plt.title(f'HG mode {int(self.nm_order[0]),int(self.nm_order[1])} intensity profile')
            plt.show()
            
    def __str__(self):
        return 'EM mode parameters are {}'.format(self.param_dict)
       
if __name__ == '__main__':
    new_mode = Modes('mode1_config.txt')
    print(new_mode)
    new_mode.plot_field('Intensity')

  