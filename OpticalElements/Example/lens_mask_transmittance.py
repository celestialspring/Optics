# -*- coding: utf-8 -*-
"""

An example with lens 
@author: SM
"""

import sys
from pathlib import Path
src_path = str(Path.cwd().parent.parent / "BeamMode_propagation")
src2_path = str(Path.cwd().parent )
if src_path not in sys.path:
    sys.path.append(src_path)
if src2_path not in sys.path:
    sys.path.append(src2_path)
import mode_class
import opticalelements
import propagation_methods
import numpy as np

#Initialize field
new_mode = mode_class.Modes('config_mask_lens.txt')
print(new_mode)
grid_size = new_mode.grid_size[0]
wavelength = new_mode.wavelength
field_init = new_mode.field_input
new_mode.plot_field('Intensity')
x,y = new_mode.generate_object_plane()

el_obj = opticalelements.OpticalElements()

#amplitude transmittance of image mask
grid_arr = x
img_arr, bin_arr = el_obj.image_amplitude_mask('../cat.jpg',grid_arr,grid_size,3E-3)

#field after mask
field_after_mask = new_mode.field_transmittance(field_init, bin_arr)

#propagate the field 
field_propagated = propagation_methods.PropagationMethods('config_mask_lens.txt', field_after_mask, x, y)

#focal length of lens
f = 20E-3

#propagate field 2f distance towards lens
field_after_2f = field_propagated.rayleigh_sommerfeld_DI(2*f)

#calculate lens phase change or transmittance
lens_transmittance = el_obj.paraxlensphase(wavelength,f,x, y)

#Apply lens transmittance to the field
field_after_lens = new_mode.field_transmittance(field_after_2f, lens_transmittance)

#Propagate the field to 2f after lens
field_after_lens_at_2f = propagation_methods.PropagationMethods('config_mask_lens.txt', field_after_lens, x, y)

field_out =field_after_lens_at_2f.rayleigh_sommerfeld_DI(2*f)

#inverted image
field_after_lens_at_2f.plot_field('Amplitude')
field_after_lens_at_2f.plot_field('Intensity of input field', field_out)



