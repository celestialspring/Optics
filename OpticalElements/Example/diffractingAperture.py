# -*- coding: utf-8 -*-
"""

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

#Initialize field
new_mode = mode_class.Modes('mode1_config.txt')
print(new_mode)
field_init = new_mode.field_input
new_mode.plot_field('Intensity')
x,y = new_mode.generate_object_plane()

#Specifying rectangular aperture
slit_widthx = 1E-4
slit_widthy = 0.2E-4
rect_obj = opticalelements.Opticalelements()
#Circular slit parameters
slit_d = 1E-4
circ_obj = opticalelements.Opticalelements()


#specify transmittance function either rectangular or circular (pass circ_t array to transmittance)
xt, yt = rect_obj.rectslit(x,y, slit_widthx, slit_widthy)
circ_t = circ_obj.circslit(x, y, slit_d )

transmittance = xt*yt

#transmitted field calculated with aperture transmittance function and initial field
field_out = new_mode.field_transmittance(field_init, transmittance)

#farfield propagation
mode_propagated = propagation_methods.PropagationMethods('mode1_config.txt', field_out, x, y)
mode_propagated.angular_spectrum()
mode_propagated.plot_field('Intensity norm xprofile')
mode_propagated.plot_field('Intensity norm')