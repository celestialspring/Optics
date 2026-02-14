# -*- coding: utf-8 -*-
"""
Created on Sat Jun 14 12:47:47 2025
A class for calculating diffraction limits based on common quantities 
@author: SM
"""
import numpy as np

class Optics_toolbox:
    def __init__(self, wavelength, wave_units):
        
        self.wavelength = wavelength
        self.wave_unit = wave_units
        self.scale_factor = self.unit_conversion(wave_units)
        
    
    def unit_conversion(self, input_units):
        """
        convert input units to a metric unit or vice versa

        Parameters
        ----------
        input_units : str 
            DESCRIPTION.
        output_units : str 
             DESCRIPTION.   
        Returns
        -------
        None.

        """
        
        metric_unit_dict = {'km':3, 'm':0, 'cm': -2, 'mm': -3, 'um': -6, 
                                  'nm':-9, 'pm':-12, 'fm':-15 }
        
        unit_in = input_units.lower()
    
        
        while True:
            if unit_in not in metric_unit_dict:
                print(f"Invalid input: expected units are {list(metric_unit_dict.keys())}")
                return None
            scale_factor  = metric_unit_dict[unit_in]
            return scale_factor
        
            
    def Abbe_limit(self, numerical_aperture):
        """
        Theoretical diffraction limit: 
        conventional (diffraction-limited) optical systems like light microscopes.
        numerical_aperture : float
           
        """
        scaling = 10**self.scale_factor
        d_limit = (self.wavelength/ (2* numerical_aperture)) * scaling 
       
        d_limit_input_unit = round(d_limit/scaling,3)
        
        return d_limit_input_unit, f"{self.wave_unit}"
    
    def Rayleigh_criterion(self, numerical_aperture):
        """
       Practical limit, for calculating Airy disk size
       numerical_aperture : float

        """
        scaling = 10**self.scale_factor
        d_limit = (1.22 * self.wavelength/ (numerical_aperture)) * scaling 
        
        d_limit_input_unit = round(d_limit/scaling,3) 
      
        return d_limit_input_unit, f"{self.wave_unit}"
    
    
    def Geometrical_NA(self, diameter, focal_length, units):
        """
        calculate NA of a simple lens
        Parameters
        ----------
        diameter : float
        Diameter of lens 
        focal_length : float
        focal_length of lens
        units: list of strings
        metric units allowed in the class Optics_toolbox
        in the order of input quantities
        Returns
        -------
        None.

        """
        if type(units) == list or type(units) == str:
            
            if type(units) == str:
                lens_NA = diameter/(2*focal_length)
                lens_NA = round(lens_NA,3)
            else: 
                
                    
                diameter_scale = self.unit_conversion(units[0])
                focale_length_scale = self.unit_conversion(units[1])
                scaling = 10**(diameter_scale-focale_length_scale)
                lens_NA = (diameter/(2*focal_length)) * scaling
                lens_NA = round(lens_NA,3)
            return lens_NA
        else:
            return None
        
       
            
    def Effective_NA(self, diameter_lens, diameter_beam, focal_length, units):
        """
        calculate effective NA of lens - beam system given their sizes

        Parameters
        ----------
        diameter_lens : float
           Diameter of lens 
        diameter_beam : float
           Diameter of Beam
        focal_length : float
           focal length of lens
        units : str or list of strings
            DESCRIPTION.

        Returns
        -------
        effective_NA : float
            DESCRIPTION.

        """
        
        if type(units) == list or type(units) == str:
            
            if type(units) == str:
                effective_NA = min(diameter_lens, diameter_beam)/(2*focal_length)
                effective_NA = round(effective_NA,3)
            else: 

                diameter_lens_scale = self.unit_conversion(units[0])
                diameter_beam_scale = self.unit_conversion(units[1])
                focale_length_scale = self.unit_conversion(units[2])
                
                if diameter_lens_scale == diameter_beam_scale:
                    effective_diameter = min(diameter_lens, 
                                               diameter_beam)
                    scaling = 10**(diameter_lens_scale -focale_length_scale)
                    effective_NA = (effective_diameter/(2*focal_length) )* scaling
                else:
                    
                    effective_diameter_scale = min(diameter_lens_scale, 
                                                   diameter_beam_scale)
                    
                    effective_diameter = min(diameter_lens*10**diameter_lens_scale, 
                                             diameter_beam*10**diameter_beam_scale)
                    
                    scaling = 10**(effective_diameter_scale -focale_length_scale)
                    effective_NA = (effective_diameter/(2*focal_length*10**focale_length_scale) )
                effective_NA = round(effective_NA,5)
            return effective_NA
        else:
            return None
        
    def Gaussian_spot_size(self, numerical_aperture):
        """
        Calculate gaussian beam spot given known NA  

        Parameters
        ----------
        numerical_aperture : float
            DESCRIPTION.

        Returns
        -------
        d_limit_input_unit : TYPE
            DESCRIPTION.
        str
            DESCRIPTION.

        """
        
        scaling = 10**self.scale_factor
        d_limit = ((2*self.wavelength) / (np.pi*numerical_aperture))* scaling
        d_limit_input_unit = round(d_limit/scaling,3)
        
        return d_limit_input_unit, f"{self.wave_unit}"
    
    def Effective_beam_spot_size(self, diameter_lens, diameter_beam, focal_length, units):
        """
        
        Calculates simple effective spot sizes given lens, beam specifications
        
        Parameters
        ----------
        diameter_lens : float
           Diameter of lens 
        diameter_beam : float
           Diameter of Beam
        focal_length : float
           focal length of lens
        units : str or list of strings
            DESCRIPTION.

        Returns
        -------
       output_dict: dictionary
       Gives different spot sizes with the available formulas
        """
        
        effective_NA = self.Effective_NA(diameter_lens, diameter_beam, focal_length, units)
        
        spot_with_abbe_limit = self.Abbe_limit(effective_NA)
        spot_with_rayleigh_limit = self.Rayleigh_criterion(effective_NA)
        spot_with_gaussian_beam = self.Gaussian_spot_size(effective_NA)

        output_dict = {'with_Abbe_limit': spot_with_abbe_limit,
        'with_Rayleigh_limit':spot_with_rayleigh_limit, 
        'Gaussian_beam_spot':spot_with_gaussian_beam}
        
        return output_dict
        