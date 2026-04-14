# -*- coding: utf-8 -*-
"""


@author: SM
A GUI to calculate relevant parameters for a fabry-perot cavity.
It imports two classes.
One which calculates reflection and transmission coefficient of a dielectric stack and another that calculates FP parameters.
The FP parameters can also be calculated for arbitrary inputs. 
Here the GUI is initialized with reflection and transmission coefficients calculated via TMM class 
for a stack of layers specified in configuration file. 

"""

import tkinter as tk
import numpy as np
from Modules import transfermat_fresnel
from Modules import fabryperot

import tkinter as tk




def load_init_param():
    try:
        TMM_class = transfermat_fresnel.TransferMatFresnel('params.txt')
        layer = TMM_class.multilayer_TMM('TM')
        r_init = round(np.real(layer['r']), 4)
        t_init = round(np.real(layer['t']), 4)
        return r_init, t_init
    except Exception as e:
        print(f"Error loading TMM logic: {e}")
        return 0.0, 0.0


def setup_gui():

    global r_calc, t_calc, wavelength, cavitylength 
    global f_val, f_entry, cf_val, cf_entry, tr_val, tr_entry, root

    root = tk.Tk()
    tab_label = tk.Label(text='Fabry-Perot GUI')
    tab_label.grid()


    rtframe = tk.Frame(root, bg='grey', borderwidth=2, 
                                 relief=tk.GROOVE, padx=2, pady=2) 
    rtframe.grid(row=1, column=0,  sticky="NSEW") 

    r_label = tk.Label(rtframe,text='Reflection coefficient (r):') 
    r_label.grid(row=0, column=0, sticky='NSEW') 
    r_calc = tk.DoubleVar()
    r_entryref = tk.Entry(rtframe,  textvariable= r_calc, state='normal', width=10) 
    r_entryref.grid(row=0, column=1,sticky='NSEW') 

    t_label = tk.Label(rtframe,text='Transmission coefficient (t):') 
    t_label.grid(row=0, column=2, sticky='NSEW') 
    t_calc = tk.DoubleVar()
    t_entryref = tk.Entry(rtframe,  textvariable= t_calc, state='normal', width=10) 
    t_entryref.grid(row=0, column=3,sticky='NSEW') 


    paramin_frame = tk.Frame(root, bg='blue', borderwidth=2, 
                                 relief=tk.GROOVE, padx=2, pady=2) 
    paramin_frame.grid(row=2, column=0, sticky='NSEW') 

    wave_label = tk.Label(paramin_frame, text='Wavelength (m):') 
    wave_label.grid(row=0, column=0, sticky='NSEW') 
    wavelength = tk.DoubleVar() 
    wave_entry = tk.Entry(paramin_frame, textvariable = wavelength, state='normal', width=10)  
    wave_entry.grid(row=0, column=1, sticky='NSEW') 

    cavitylength = tk.DoubleVar() 
    cav_label = tk.Label(paramin_frame, text='Cavity length (m):') 
    cav_label.grid(row=0, column=2, sticky='NSEW') 
    cav_entry = tk.Entry(paramin_frame, textvariable = cavitylength, state='normal', width=10)  
    cav_entry.grid(row=0, column=3, sticky='NSEW') 


    result_frame = tk.Frame(root, bg='green', borderwidth=2, 
                                 relief=tk.GROOVE, padx=2, pady=4) 
    result_frame.grid(row=3, column=0, sticky='NSEW') 

    f_label = tk.Label(result_frame, text='Finesse:') 
    f_label.grid(row=0, column=0, sticky='NSEW') 
    f_val = tk.StringVar() 
    f_entry = tk.Entry(result_frame, textvariable=f_val, state='readonly', width=10) 
    f_entry.grid(row=0, column=1, sticky='NSEW') 

    cf_label = tk.Label(result_frame, text='Coeff of Finesse:') 
    cf_label.grid(row=0, column=2, sticky='NSEW') 
    cf_val = tk.StringVar() 
    cf_entry = tk.Entry(result_frame, textvariable=cf_val, state='readonly', width=10) 
    cf_entry.grid(row=0, column=3, sticky='NSEW') 

    tr_label = tk.Label(result_frame, text='Transmittance:') 
    tr_label.grid(row=0, column=4, sticky='NSEW') 
    tr_val = tk.StringVar() 
    tr_entry = tk.Entry(result_frame, textvariable=tr_val, state='readonly', width=10) 
    tr_entry.grid(row=0, column=5, sticky='NSEW') 


    calculate_frame = tk.Frame(root, bg='red', borderwidth=2, 
                                 relief=tk.GROOVE, padx=2, pady=2) 
    calculate_frame.grid(row=1, column=1, rowspan=3, sticky='NSEW') 

    calculate_button = tk.Button(calculate_frame, text='Compute', width=9, height=3, bg='gray', command=compute_results)
    calculate_button.grid(row=0, column=0, sticky='NSEW')

    return root

def compute_results():
    try:
        r = float(r_calc.get())
        lambda_ = wavelength.get()
        cavlen = cavitylength.get() 
        
        fbclass = fabryperot.FabryPerot(lambda_, r, cavlen)
        F = fbclass.coeff_Of_finesse()
 
        cf_entry.config(state='normal')
        cf_val.set(round(F, 4))
        cf_entry.config(state='readonly')
        
        Fi = fbclass.finesse()
        
        f_entry.config(state='normal')
        f_val.set(round(Fi,4))
        f_entry.config(state='readonly')
        
        Transmittance = fbclass.transmittance(F)
        tr_entry.config(state='normal')
        tr_val.set(round(Transmittance,4))
        tr_entry.config(state='readonly')
        
    except Exception as e:
        print(f"Calculation Error: {e}")


if __name__ == '__main__':

    main_window = setup_gui()
    r_val_start, t_val_start = load_init_param()
    
    r_calc.set(r_val_start)
    t_calc.set(t_val_start)
    
    main_window.mainloop()


