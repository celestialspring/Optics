# -*- coding: utf-8 -*-
"""
Gaussian Beam Propagation (free space)

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
    
    def gaussian_translation_mat(self,dist):
        
        T = np.array([[1, dist],
                      [0, 1]])
        
        return T
    
    def gaussian_refraction_mat(self, focallength):
        
        power = 1/focallength
        R = np.array([[1, 0],
                      [-power, 1]])
        
        return R
    
    def gaussian_single_lens_transform(self, obj_space_dist, focallength, img_space_dist):
        self.lens_position(obj_space_dist, focallength, img_space_dist)
        M = self.gaussian_translation_mat(img_space_dist) \
            @ self.gaussian_refraction_mat(focallength) \
            @ self.gaussian_translation_mat(obj_space_dist)
            
        A = M[0,0]
        B = M[0,1]
        C = M[1,0]
        D = M[1,1] 
        
        return A, B, C, D
    
    def lens_position(self, dist_from_waist, focal_length, imaging_dist):
        
        self.lens_dist = dist_from_waist
        self.focal_length = focal_length
        self.imaging_dist = imaging_dist
        
    def gaussian_ABCD(self,q,A,B,C,D):
        
        q_k = (A*q+B)/(C*q+D)
        self.q = q_k
        return self.q
    
    def param_q(self, q):
        
        z_dist_to_waist = np.real(q)
        self.z_r = np.imag(q)
        return z_dist_to_waist, self.z_r
    
    def update_params(self, waist=0, div_angle=0, z_r=0):
        self.waist = waist
        self.div_angle = div_angle
        self.z_r = z_r
        
    def define_q(self, z_dist_to_waist, z_r):
        
        self.q = z_dist_to_waist + 1j*z_r
        return self.q
    
    def simulation_plane(self):
        
        lx = 10 * self.waist
        lz = 3* self.z_r
        nz = 500 
        nx = 200
        x =  np.linspace(0, lx, nx)-lx/2
        z = np.linspace(0, lz, nz) 
        
        
        Z, X = np.meshgrid(z,x)
        
        return X,Z
    
    def segmented_simulation_plane(self, object_distance, image_distance):
        N_obj = object_distance/self.wavelength
        N_image = image_distance/self.wavelength
        lx = 300* self.wavelength
        lz =  1.5*(N_obj+N_image)* self.wavelength
        
       
        nz_total = 1000 
        z_global = np.linspace(0, lz, nz_total)
        
       
        z1 = z_global[z_global <= object_distance]
        z2 = z_global[z_global > object_distance]
        
        nx = 200
        x = np.linspace(0, lx, nx) - lx/2
        
        
        Z1, X1 = np.meshgrid(z1, x)
        Z2, X2 = np.meshgrid(z2, x)
        return x,X1,X2,Z1,Z2
    
        
    def plot(self, plottype:str, plotvar,X,Z):
        
       
        if plottype == 'free space propagation':
            norm_z = Z/self.z_r
            norm_x = X/self.waist
            Intensity = plotvar
            plt.figure()
            plt.imshow(Intensity, cmap ='inferno', extent=[np.min(norm_z), np.max(norm_z), np.min(norm_x), np.max(norm_x)], aspect='auto')
            plt.colorbar(label='Intensity')
            plt.xlabel('$z/z_R$')
            plt.ylabel('$x/w_0$')
            plt.title("Gaussian beam propagation")
    
            info_text = f"$z_r$ = {round(self.z_r,2)} um, $w_0$ = {round(self.waist,4)} um"
            plt.text(0.25, 0.95, info_text, transform=plt.gca().transAxes, 
                     fontsize=10, verticalalignment='top', color='white')
    
        if plottype == 'single lens propagation':
            norm_z = Z
            norm_x = X
            Intensity_lens = plotvar
            plt.figure()
            plt.imshow(Intensity_lens, cmap ='inferno', extent=[np.min(norm_z), np.max(norm_z), np.min(norm_x), np.max(norm_x)], aspect='auto')
            plt.colorbar(label='Intensity')
            plt.xlabel('$z/\lambda$')
            plt.ylabel('$x/\lambda$')
            plt.title("Gaussian beam propagation with single lens")
    
            info_text = f"$z_r$ = {round(self.z_r,2)} mm, $w_0$ = {round(self.waist,4)} mm"
            
            plt.text(0.02, 0.95, info_text, transform=plt.gca().transAxes, 
                     fontsize=10, verticalalignment='top', color='white')
            plt.axvline(self.lens_dist/self.wavelength+.1, color='white', linestyle='-', linewidth=1.5)
            plt.text(self.lens_dist/self.wavelength+0.1, np.min(norm_x), 'Lens', color='white', fontsize=10, 
                 verticalalignment='bottom')
            plt.axvline((self.lens_dist+self.focal_length)/self.wavelength, color='green', linestyle='--', linewidth=1.5)
            plt.text((self.lens_dist+self.focal_length)/self.wavelength, -100, 'Focal plane', color='green', fontsize=10, 
                 verticalalignment='bottom')
            if self.focal_length < self.lens_dist:
                plt.axvline((self.lens_dist-self.focal_length)/self.wavelength, color='green', linestyle='--', linewidth=1.5)
            plt.axvline((self.lens_dist+self.imaging_dist)/self.wavelength, color='red', linestyle='--', linewidth=1.5)
            plt.text((self.lens_dist+self.imaging_dist)/self.wavelength, 100, 'Image plane', color='red', fontsize=10,
                     verticalalignment='bottom')
            
if __name__ == '__main__':
    beam = GaussianBeamPropagation(632E-9, 0.5E-3)
    param = beam.params
    x, z = beam.simulation_plane()
    intensity = beam.gaussian_field(x,z)
    print(intensity.shape)
    beam.plot('free space propagation', intensity, x,z)
  