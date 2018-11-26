from __future__ import print_function
from functools import partial
import numpy as np
from configparser import ConfigParser
import iminuit
import matplotlib.pyplot as plt

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
    
def minimize(likelihood, model, wave, flux, ivar, **init_pars):
    like = partial(likelihood, line_model=model, wave=wave, flux=flux, ivar=ivar)

    par_names = [p for p in model.__code__.co_varnames[:model.__code__.co_argcount]]
    #print(par_names)

    m = iminuit.Minuit(like,
        forced_parameters = par_names,
        errordef = 1, pedantic = False,
        **init_pars)
    
    fmin  = m.migrad()
    return m,fmin

def double_minimize(likelihood1, likelihood2, model, wave, flux, ivar, **init_pars1):
    
    par_names = [p for p in model.__code__.co_varnames[:model.__code__.co_argcount]]
    
    m1, fmin1 = minimize(likelihood1, model, wave, flux, ivar, **init_pars1)
    
    init_pars2 = {x:y for x,y in zip(par_names, m1.values.values())}

    m2, fmin2 = minimize(likelihood2, model, wave, flux, ivar, **init_pars2)
    return m1, fmin1, m2, fmin2

def plot_fit(wave, model, m, color, lab):
    plt.plot(wave, model(*[m.values[p] for p in m.parameters], wave=wave), color, lw=2, alpha = .8, label = lab)
    
def get_chi(like, model, m, wave, flux, ivar):
    return str(like(*[m.values[p] for p in m.parameters], line_model = model, wave=wave, flux = flux, ivar = ivar)) +' / ' + str((len(wave) - len(m.parameters)))

    
