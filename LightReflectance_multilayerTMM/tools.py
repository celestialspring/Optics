# -*- coding: utf-8 -*-
"""

@author: SM

A simple function to convert txt file of a given format into python dictionary
"""
import os

def read_configtxt(config_file:str):
    current_wd = os.getcwd()
    file_path = current_wd + f'\{config_file}' 
    if not os.path.exists(file_path):
        raise Exception('The config_file does not exists.')
    else:
        file = open(config_file, encoding="utf8") #utf8 for non-standard characters 
        lines = file.readlines()
        file.close()
        
    config_dict = {}    
    
    for line in lines:
        #clean variables and add to dictionary
        line_list=line.strip().split(':')
        key = line_list[0].strip(" '' ")
        value_list = line_list[1].split('#')
        val = value_list[0].strip()
        if '(' and ')' in val:
            clean_val = val.strip("()")
            x, y = clean_val.split(",")
            val = (float(x), float(y))
        else:
            try:
                val = float(val)
            except:
                pass
        config_dict[key] = val
   
    
    return config_dict

if __name__ =='__main__':
   read_configtxt('mode1_config.txt')
   