# -*- coding: utf-8 -*-
"""
This script simulates propagation of two Hermite-Gaussian modes specified in two config files.
It calls methods from two classes: mode_class and propagation_methods

@author: SM
"""
import mode_class
import propagation_methods

'Input files for two modes'
mode1_config_file = 'mode1_config.txt'
mode2_config_file = 'mode2_config.txt'

'Initial field and computational grid for the first mode'
mode_1 = mode_class.Modes(mode1_config_file)
X1, Y1 = mode_1.generate_object_plane()
mode1_init_field = mode_1.field_input

'Input field and grid fro the second mode'
mode_2 = mode_class.Modes(mode2_config_file)
X2, Y2 = mode_2.generate_object_plane()
mode2_init_field = mode_2.field_input

'Propagate the first mode and get the output field after FFT-DI'
propagate_mode1 = propagation_methods.PropagationMethods(mode1_config_file,
                                                          mode1_init_field, X1, Y1)
mode1_field_final = propagate_mode1.rayleigh_sommerfeld_DI()

'Propagate the second mode and get its field after FFT-DI'
propagate_mode2 = propagation_methods.PropagationMethods(mode2_config_file,
                                                          mode2_init_field, X2, Y2)
mode2_field_final = propagate_mode2.rayleigh_sommerfeld_DI()

'Calculate the overlap between the two output fields for the same grid size'
overlap = propagation_methods.two_mode_overlap(mode1_field_final, mode2_field_final, X1, Y1)
print(f'The overlap between two modes is {overlap}')

'Profiles and 2D maps of modes'
mode_1.plot_field('Intensity xprofile')
propagate_mode1.plot_field('Intensity xprofile')
propagate_mode1.plot_field('Intensity')
mode_2.plot_field('Intensity xprofile')
propagate_mode2.plot_field('Intensity xprofile')
propagate_mode2.plot_field('Intensity')