# -*- coding: utf-8 -*-
"""
Created on Sun Dec 28 14:42:22 2025

@author: SandyM
"""

# HG mode propagation

This project attempts to simulate free space propagation of Hermite-Gaussian mode 
using fast fourier transform (FFT) based direct integration method (DI) 
of Rayleigh-Sommerfeld diffraction integral. 
It was created as part of Cailabs exercise. 

## Objective
- To compute transverse profiles of a set of 2 Hermite-Gaussian modes of given orders
  (m1, n1) and (m2, n2) after numerical propagation over a given distance in free space.
- Return a scalar value corresponding to the overlap between two modes after propagation

## Scripts
- `simulate_modes.py`: Run this script to generate outputs for inputs given in .txt files 
- `propagation_methods.py`: Contains a class used for implementing FFT-DI as well as a function to evaluate mode overlap
- `mode_class.py`: Contains a class used for generating input field of the specified mode and object plane grid
- `tools.py`: Contains a function that converts .txt file of a certain format into python dictionary

## General architecture 

### Python libraries 
- os, numpy, matplotlib, scipy, copy, warnings

### Inputs
An example of the inputs required in a .txt file is shown in the following:
```text
'mode_type':           Hermite-Gaussian #
'wavelength':          1000E-9 #in meters
'field_amplitude':     1 #constant amplitude of field
'waist_radius':        150E-6 #in meters
'n_m':		       (2,0) #n,m order of hermite polynomial upto(3,3)
'sample_points':       600 #
'grid_size':           (600E-6,600E-6) # in meters
'propagation_distance': 5E-2 #in meters
```
### Implementation 
The complete numerical implementation follows the structure given in the script `simulate_modes.py` 
- Indicate input .txt file
- Create mode object from the class `Mode()`in `mode_class.py` with the input .txt file
- Pass the mode field input and X,Y arrays defining the grid or object plane 
  to `PropgationMethods()` class in `propagation_methods.py`
- The field after propagation is returned by calling `rayleigh_sommerfeld_DI()` function on the object of `PropgationMethods()`
- The above procedure is repeated for a second mode 
- Finally, the overlap of the two modes is calculated by passing mode fields and X, Y arrays that defines a grid 
  into `two_mode_overlap()` function in `propagation_methods.py`
  
### Outputs
In the script `simulate_modes.py`, the outputs are:
- Transverse profile plots corresponding to modes before and after propagation as well as
  a 2D image of the mode after propagation are plotted
- Scalar value of overlap integral

### How to run
From a terminal, once in the directory of the script, run:
```python

simulate_modes.py
``` 
Note: Similarly all .py files can be run as a standalone script. 

### Limitations and edge cases 
- The .txt file format is fixed. It cannot handle multiple values, presence of two colons etc.
- The .txt files must be in the same directory as the main scripts.
- HG modes only available upto order (3,3).
- If propagation distance is zero, it returns initial field.
- If z too close, must increase points (rapid oscillations in near field). Increasing window smoothens out    	everything
- No attempt at optimization is made. Test window 50-800 points. 
- Sampling points must be the same between two modes for calculating their overlap with `two_mode_overlap()`.
- The overlap calculation considers the same computational grid area for both the modes.
   
