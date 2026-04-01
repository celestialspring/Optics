# -*- coding: utf-8 -*-
"""

@author: SM
"""

# Tranfer matrix approach and Fresnel equations 

The script `transfermat_fresnel.py` contains a class `TransferMatFresnel()` that has methods to calculate transmission and reflection coefficients. 


## Objective
- Compute r, t, R, T values for different layer stacks or interfaces using TMM or Fresnel equations.

## Scripts
- `transfermat_fresnel.py`: Run this script to generate outputs for inputs given in .txt file
- `tools.py`: Contains a function that converts .txt file of a certain format into python dictionary

## General architecture 

### Python libraries 
- numpy

### Inputs
An example of the inputs required in a .txt file is shown in the following:
```text
'wavelength':          600E-9 #in meters
'layers':	       3 #int number of layers excluding substrate, first layer by default is air (n0=1)
'incident_angle':      45 #in degree
'nsubstrate':          1.5 # float, n of substrate
'n1':                  2.4     #float
't1':                  50E-9   #layer1 thickness in meters
'n2':                  1.45    #float
't2':                  100E-9  #layer2 thickness in meters
'n3':                  2.4     #float
't3':                  50E-9   #layer3 thickness in meters
```
### Functions
- loadparam(self, configfile): takes in .txt file as a string and initializes variables 
- singlelayer_TMM(self, polarisation: str): evaluate coefficients of a single layer, polarisation can be TE or TM
- fresnel_interface(self, polarisation:str): evaluate coefficients at an interface for TM or TE polarisation 
- multilayer_TMM(self, polarisation:str): evaluate coefficients at the output of a stack of layers as given in the example above
  
### Outputs
- A dictionary with keys r, t, R, T and corresponding values

### How to run
From a terminal, once in the directory of the script, run:
```python
transfermat_fresnel.py
``` 

   
