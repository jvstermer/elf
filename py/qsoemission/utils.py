from __future__ import print_function
from functools import partial
import numpy as np
from configparser import ConfigParser
import iminuit

from . import const

def go_to_lf(wave, z):
    '''
    Converts wavelenghts from rest frame to laboratory frame

    Arguments:
        lam -- wavelength (numpy array)
        z -- redshift (float)

    Returns:
        log10(lab wavelength)
    '''

    return wave*(1+z)

def get_system_values(path, section, value):
    cp = ConfigParser()       
    cp.read(path)
    return cp.get(section,value)

def window(z, wave, flux, ivar, line_id, range_window):
    '''
    Restrict the data to a window around an emission line

    Arguments:
        z -- quasar redshift (float)
        flux -- quasar flux (nupy array)
        wave -- wavelength in laboratory frame (numpy array)
        ivar -- inverse variances (numpy array)
        line_id -- name of the line (str)
    
    Returns:
        wavelength -- wavelength restricted to the window (numpy array)
        flux -- flux restricted to the window (numpy array)
        ivar -- ivar restricted to the window (numpy array)
    '''

    line = const.emission_lines[line_id]
    mask = (wave >= go_to_lf(line-range_window/2, z)) & (wave <= go_to_lf(line+range_window/2, z))

    return wave[mask], flux[mask], ivar[mask]
    
def minimization(likelihood, model, wave, flux, ivar, **init_pars):

    like = partial(likelihood, line_model=model, wave=wave, flux=flux, ivar=ivar)
    print(list(init_pars.keys()))
    m = iminuit.Minuit(like,
        forced_parameters = list(init_pars.keys()),
        errordef = 1,
        **init_pars)

    fmin  = m.migrad()

    return m,fmin
