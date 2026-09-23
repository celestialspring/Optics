# -*- coding: utf-8 -*-
"""


@author: SM
"""

obj_space_dist = 100 #waist at f distance in object space
focallength = 50
img = 1/focallength - 1/obj_space_dist
img_space_dist = 1/img
wavelength = 0.000850
beam2 = GaussianBeamPropagation(wavelength, 0.100)
params = beam2.params
x,x1,x2, z1,z2 = beam2.segmented_simulation_plane(obj_space_dist,img_space_dist)
intensity1 = beam2.gaussian_field(x1,z1)
z_dist_to_waist = 0
w_0 = params['waist']
z_r = params['rayleigh_range']
q0 = beam2.define_q(z_dist_to_waist, z_r) 
print(f'{q0} q0 and selfq {beam2.q}')
  
A,B,C,D = beam2.gaussian_single_lens_transform(obj_space_dist, focallength, img_space_dist)
q1 = beam2.gaussian_ABCD(q0, A, B, C, D)
  
z_dist_to_newwaist, new_z_r = beam2.param_q(q1)
  
print(f'{z_dist_to_newwaist}:z_dist_to_newwaist, {new_z_r}: new_z_r,selfq {beam2.q}')
beam2.update_params(0,0,new_z_r)
beam2.gaussian_param(n=1)


  # intensity2 = beam2.gaussian_field(x2,z2-obj_space_dist-img_space_dist+z_dist_to_newwaist)
  # tot_intensity = np.hstack((intensity1,intensity2))
  # z = np.hstack((z1,z2))/wavelength

  # beam2.plot('single lens propagation', tot_intensity,x/wavelength,z)