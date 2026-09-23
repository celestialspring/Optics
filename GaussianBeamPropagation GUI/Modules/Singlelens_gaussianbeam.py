# -*- coding: utf-8 -*-
"""


@author: SM
"""
import numpy as np
import GaussianBeam

#units in mm
#Distance in the object space to lens
obj_space_dist = 100 
#f lens
focallength = 50
#imaging plane
img = 1/focallength - 1/obj_space_dist
img_space_dist = 1/img

#units in mm
wavelength = 0.000850
waist = 0.100
beam = GaussianBeam.GaussianBeamPropagation(wavelength, waist)
params = beam.params
z_dist_to_waist = 0
w_0 = params['waist']
z_r = params['rayleigh_range']
q0 = beam.define_q(z_dist_to_waist, z_r) #q parameter
print(f'q0 is {beam.q}')

#Define simulation plane based on obj and img distances  
x,x1,x2, z1,z2 = beam.segmented_simulation_plane(obj_space_dist,img_space_dist)
intensity1 = beam.gaussian_field(x1,z1)

#Propagate the beam using ABCD matrix method
A,B,C,D = beam.gaussian_single_lens_transform(obj_space_dist, focallength, img_space_dist)
q1 = beam.gaussian_ABCD(q0, A, B, C, D)
  
#new q parameter after lens
waist_to_currentdist, new_z_r = beam.param_q(q1)
print(f'{waist_to_currentdist}:waist_to_currentdistance, {new_z_r}: new_z_r')

#update object attributes associated to the new beam parameters
beam.update_params(0,0,new_z_r)
beam.gaussian_param(n=1)

#the simulation plane is segmented. To the right of the lens, the maximum of intensity
#is recentered to beam waist. The image space starts after object space:z-obj_space
#The new waist appears at x distance from image plane: img_plane-Re(q1)
  
z2_recenter = z2 - obj_space_dist-img_space_dist+waist_to_currentdist
intensity2 = beam.gaussian_field(x2,z2_recenter)

#concatenate the two simulation planes
tot_intensity = np.hstack((intensity1,intensity2))
#normalize with wavelength
Z = np.hstack((z1,z2))/wavelength
X =  x/wavelength
beam.plot('single lens propagation', tot_intensity,X,Z)