# -*- coding: utf-8 -*-
"""


@author: SM
"""
from PIL import Image 
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
    
    def paraxlensphase(self, wavelength, focallength, x_grid, y_grid,
                       lens_index=1.51, lens_diameter=10E-3, thickness=9E-3):
        
        k = (2*np.pi)/wavelength
        if focallength == 0:
            phase = np.exp(1j*k*lens_index*thickness)
        else:
            phase = np.exp(1j*k*lens_index*thickness)*np.exp(-1j*np.pi*(x_grid**2+y_grid**2)/(wavelength*focallength)) 
        
        return phase
    
    def image_amplitude_mask(self, filename:str, grid_arr, grid_size, mask_size):
        """
        Amplitude mask based on an image

        Parameters
        ----------
        filename : str
            give the filename and type of the image
        grid_arr : TYPE
            simulation grid array
        grid_size : TYPE
            simulation_window in physical dimension
        mask_size : TYPE
            Targeted size of the mask

        Returns
        -------
        img_array : TYPE
            resized array of the mask 
        binary_arr : TYPE
            mask inside wider simulation field array

        """
        N= grid_arr.shape[0]
        dx = grid_size/N
        pixels = int(mask_size/dx)
        
        img = Image.open(filename).convert('L')
        resize_img = img.resize((pixels, pixels))

        img_array = np.array(resize_img) / 255
        binary_arr = np.zeros_like(grid_arr, dtype=(float))

        
        start_index = (N - pixels )//2
        end_index = start_index+ pixels
        binary_arr[start_index:end_index, start_index:end_index] = img_array
        return img_array, binary_arr
    
if __name__ == '__main__':
    x = np.linspace(0,8E-3,100) -4E-3
    y = np.linspace(0,8E-3,100) -4E-3

    X,Y = np.meshgrid(x,y)
    clss = OpticalElements()
    xr, yr = clss.rectslit(X,Y, 100E-5,100E-5)
    mask = xr*yr
    plt.figure()
    plt.imshow(mask, extent=[np.min(x), np.max(x), np.min(y), np.max(y)])
    plt.colorbar()
    plt.show()
    
    img_array, amplitude_mask = clss.image_amplitude_mask('cat.jpg', X, 10E-3, 0.5E-3)
    plt.figure(figsize=(8, 4))
    plt.subplot(1, 2, 1)
    plt.title("Original Greyscale")
    plt.imshow(img_array, cmap='gray')
    
    plt.subplot(1, 2, 2)
    plt.title("Binary Amplitude Mask ")
    plt.imshow(amplitude_mask, cmap='gray')
    plt.show()