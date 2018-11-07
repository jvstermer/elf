import numpy as np
from configparser import ConfigParser

def go_to_lf(lam, z):
    return np.log10( lam * (z+1) )

def set_system_values(path, section, value):

    cp = ConfigParser()       
    cp.read(path)
    return cp.get(section,value)
