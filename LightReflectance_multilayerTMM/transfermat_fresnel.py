# -*- coding: utf-8 -*-
"""

@author: SM
"""
import tools
import numpy as np


"constants"
epsilon0 = 8.854*10**(-12) #F m-1
mu0 = 4*np.pi*10**(-7) # N A-2
 
class TransfermatFresnel():
    
    def __init__(self, configfile: str):
        self.loadparam(configfile)
    
    def loadparam(self, configfile):
        
        self.paramdict = tools.read_configtxt(configfile)
        self.wavelength = self.paramdict['wavelength']
        self.layers = int(self.paramdict['layers'])
        self.incident_angle = int(self.paramdict['incident_angle'])
        self.nsub = self.paramdict['nsubstrate']
        self.nt_tupslist = []
        if self.layers >= 1:
            if not len(self.paramdict)-4 == 2*self.layers: #each layer must have two inputs n,t
                print('Data insufficient considering the case of only a single layer')
                self.nt_tupslist.append((self.paramdict['n1'],self.paramdict['t1']))
            else:
                for N in range(self.layers):
                    n_num = 'n'+str(N+1)
                    t_num = 't'+str(N+1)
                    self.nt_tupslist.append((self.paramdict[n_num],self.paramdict[t_num]))  
                    
                
    def singlelayer_TMM(self, polarisation: str):
        
        n0 = 1
        n1, t1 = self.nt_tupslist[0]
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
            gamma1 = (n1*np.sqrt(epsilon0*mu0))/np.cos(theta_radts)  
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
        
        #transmissitance and reflectance
        R= round(abs(r)**2,4)
        T= round(1- R,4)
        return {'r': r, 't':t, 'R': R, 'T': T}
    
if __name__ =='__main__':
    layer1 = TransfermatFresnel('params.txt')
    print(layer1.singlelayer_TMM('TE'))
        
        
        
        
