# -*- coding: utf-8 -*-
"""

@author: SM
"""
import tools
import numpy as np


"constants"
epsilon0 = 8.854*10**(-12) #F m-1
mu0 = 4*np.pi*10**(-7) # N A-2
 
class TransferMatFresnel():
    
    def __init__(self, configfile: str):
        self.loadparam(configfile)
    
    def loadparam(self, configfile):
        
        self.paramdict = tools.read_configtxt(configfile)
        self.wavelength = self.paramdict['wavelength']
        self.layers = int(self.paramdict['layers'])
        self.incident_angle = int(self.paramdict['incident_angle'])
        self.nsub = self.paramdict['nsubstrate']
        self.nt_tupslist = []
       
        if self.layers == 2:
            
                for N in range(self.layers):
                    n_num = 'n'+str(N)
                    t_num = 't'+str(N)
                    try:
                        self.nt_tupslist.append((complex(self.paramdict[n_num]),self.paramdict[t_num]))
                    except:
                        if n_num == 'n0':
                            print(' considering n0 as air')
                            n0 = 1
                            t0 = 0
                            self.nt_tupslist.append((complex(n0),t0))
                            continue
        elif self.layers > 2:  
    
            if not len(self.paramdict)-4 == 2*self.layers: #each layer must have two inputs n,t
                print('Data insufficient considering the case of 2 layers with air(n0)-n1')
                self.nt_tupslist.append((complex(self.paramdict['n1']),self.paramdict['t1']))
            else:
                for N in range(self.layers):
                    n_num = 'n'+str(N)
                    t_num = 't'+str(N)
                    try:
                        self.nt_tupslist.append((complex(self.paramdict[n_num]),self.paramdict[t_num]))
                    except:
                        if n_num == 'n0':
                            print(' considering n0 as air')
                            n0 = 1
                            t0 = 0
                            self.nt_tupslist.append((complex(n0),t0))
                            continue
        else:
            raise Exception('Number of layers must be at least 2')        
                    
                
    def singlelayer_TMM(self, polarisation: str):
        '''
        Transfer matrix method to calculate the R and T of a single thin layer comparable
        to wavelength
        Parameters
        ----------
        polarisation : str
            TE or TM
    

        Returns
        -------
        dict
            DESCRIPTION.

        '''
        
        n0, t0 = self.nt_tupslist[0]
        n1, t1 = self.nt_tupslist[1]
        k0 = (2*np.pi)/self.wavelength
        
        
        theta_radt1 = np.arcsin((n0*np.sin(np.radians(self.incident_angle)))/n1)
        theta_radts = np.arcsin((n1*np.sin(theta_radt1))/self.nsub)
        delta =  n1*k0*t1*np.cos(theta_radt1)
    
        if polarisation == 'TE':
            gamma0 = n0*np.sqrt(epsilon0*mu0)*np.cos(np.radians(self.incident_angle))
            gamma1 = n1*np.sqrt(epsilon0*mu0)*np.cos(theta_radt1)
            gammas = self.nsub * np.sqrt(epsilon0*mu0)*np.cos(theta_radts)
        else:
            gamma0 = n0*np.sqrt(epsilon0*mu0)/np.cos(np.radians(self.incident_angle))
            gamma1 = (n1*np.sqrt(epsilon0*mu0))/np.cos(theta_radt1)  
            gammas = self.nsub * np.sqrt(epsilon0*mu0)/np.cos(theta_radts)
        
        m11 = np.cos(delta)
        m12 = (1j*np.sin(delta))/gamma1
        m21 = 1j*gamma1*np.sin(delta)
        m22 = np.cos(delta)
        
        t_num = 2*gamma0
        denom = gamma0*m11 + gamma0*gammas*m12+ m21+gammas*m22
        r_num = gamma0*m11 + gamma0*gammas*m12- m21-gammas*m22
        
        #transmission and reflection coefficients
        t= t_num/denom
        r = r_num/denom
        
        #transmissitance and reflectance #coherent case or t~lambda
        R= round(abs(r)**2,4)
        T= round(1- R,4)
   
        return {'r': r, 't':t, 'R': R, 'T': T}
    
    def fresnel_interface(self, polarisation:str):
        '''
        Fresnel calculation to estimate T and R at the boundary
        Parameters
        ----------
        polarisation : str
            TE or TM

        Returns
        -------
        dict
            DESCRIPTION.

        '''
     
        n0, t0 = self.nt_tupslist[0]
        n1, t1 = self.nt_tupslist[1]
        k0 = (2*np.pi)/self.wavelength
 
        theta_radt1 = np.arcsin((n0*np.sin(np.radians(self.incident_angle)))/n1)

      
        if polarisation == 'TE':
           r_num = n0*np.cos(np.radians(self.incident_angle))-n1*np.cos(theta_radt1)
           r_denom =  n0*np.cos(np.radians(self.incident_angle))+n1*np.cos(theta_radt1)
           
           r = r_num/r_denom
           t = 1 + r
        else:
           r_num = -n1*np.cos(np.radians(self.incident_angle))+n0*np.cos(theta_radt1)
           r_denom =  n1*np.cos(np.radians(self.incident_angle))+n0*np.cos(theta_radt1)
           
           r = r_num/r_denom
           t = (n0*(1 - r))/n1
    
        R= round(abs(r)**2,4)
        T= round(1- R,4)
        
        return {'r': r, 't':t, 'R': R, 'T': T}

    def multilayer_TMM(self, polarisation:str):
        
        n0, t0 = self.nt_tupslist[0]
        k0 = (2*np.pi)/self.wavelength
        S = n0*np.sin(np.radians(self.incident_angle))
        M_total = np.array([[1.0 + 0j, 0.0 + 0j], 
                        [0.0 + 0j, 1.0 + 0j]])
        
        if polarisation == 'TE':
            gamma0 = n0*np.sqrt(epsilon0*mu0)*np.cos(np.radians(self.incident_angle))
            gammas = np.sqrt(epsilon0*mu0)* np.sqrt(self.nsub**2-S**2)
        else:
            gamma0 = n0*np.sqrt(epsilon0*mu0)/np.cos(np.radians(self.incident_angle))
            gammas = self.nsub**2 * np.sqrt(epsilon0*mu0)/np.sqrt(self.nsub**2-S**2)
        
        if self.layers >= 2:
            if self.layers == len(self.nt_tupslist):
                for i in range(len(self.nt_tupslist)):
                    if i == 0:
                        continue
                    else:
                        n,t = self.nt_tupslist[i]
                        if polarisation == 'TE':
                            gamma_j = np.sqrt(epsilon0*mu0)*np.sqrt(n**2 - S**2)
                        else: 
                            gamma_j = n**2*np.sqrt(epsilon0*mu0) / (np.sqrt(n**2 - S**2))
                        delta_j = k0*t*np.sqrt(n**2-S**2)
                        M11 = np.cos(delta_j)
                        M12 = 1j*np.sin(delta_j)/(gamma_j)
                        M21 = 1j*gamma_j*np.sin(delta_j)
                        M22 = np.cos(delta_j)
                        M_j = np.array([[M11,M12],
                                        [M21, M22]])
                        M_total = M_total @ M_j
                                        
                M11_tot = M_total[0, 0]
                M12_tot = M_total[0, 1]
                M21_tot = M_total[1, 0]
                M22_tot = M_total[1, 1]
    
               
                t_num = 2 * gamma0
                denom = (gamma0 * M11_tot) + (gamma0 * gammas * M12_tot) + M21_tot + (gammas * M22_tot)
                r_num = (gamma0 * M11_tot) + (gamma0 * gammas * M12_tot) - M21_tot - (gammas * M22_tot)
                
                #transmission and reflection coefficients
                t= t_num/denom
                r = r_num/denom
                
                #transmissitance and reflectance #coherent case or t~lambda
                R= round(abs(r)**2,4)
                T= round(1- R,4)
   
        return {'r': r, 't':t, 'R': R, 'T': T}
        
   

if __name__ =='__main__':
    layer1 = TransferMatFresnel('params.txt')
    # print(layer1.singlelayer_TMM('TM'))
    # print(layer1.fresnel_interface('TM'))
    print(layer1.multilayer_TMM('TM'))

        
        
        
        
