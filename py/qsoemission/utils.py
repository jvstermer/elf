from __future__ import print_function
from functools import partial
import numpy as np
from configparser import ConfigParser
import iminuit
import matplotlib.pyplot as plt
from py.qsoemission import likelihood, line_models

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

def pick_method(method):
    if method != "chi_squared": 
        like_0 = likelihood.chi_squared
    if method == "chi_squared":
        like_0 = None
    return like_0

def get_pars(model):
    return [p for p in model.__code__.co_varnames[:model.__code__.co_argcount]]

def get_parnames(model_label, system, model):
    
    if model_label == 'spl':
        nodes = int(get_system_values(system, 'line model', 'nodes'))
        pars_model = ['par_{}'.format(i) for i in range(nodes)]
        noise = None
        pars_noise = []
        noise_label = None
    else:
        pars_model = get_pars(model)
        noise_label = get_system_values(system, 'bkg model', 'bkg')
        noise = getattr(line_models, noise_label)
        pars_noise = get_pars(noise)
        noise = line_models.line_model(noise, pars_noise)
        
    model = line_models.line_model(model, pars_model)
    return model, noise, noise_label, pars_model, pars_noise
        
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
        
def minimize(likelihood, line, model, wave, flux, ivar, noise=None, x = None, **init_pars):

    like = partial(likelihood, line = line, wave = wave, flux = flux, ivar = ivar, x = x, model = model, noise = noise)
    
    par_names = line.parnames
    
    m = iminuit.Minuit(like, 
                       forced_parameters = par_names, 
                       errordef = 1, pedantic = False, 
                       **init_pars)
    fmin  = m.migrad()
    return m, fmin

def double_minimize(likelihood1, line, model,  wave, flux, ivar, noise=None, x = None, likelihood2 = None, **init_pars):
    if likelihood2 == None:
        m, fmin = minimize(likelihood1, line, model, wave, flux, ivar, noise, x, **init_pars)
    
    else:
        m, fmin = minimize(likelihood2, line, model, wave, flux, ivar, noise, x, **init_pars)
        par_names = line.parnames
        init_pars2 = {x:y for x,y in zip(par_names, m.values.values())}
        
        m, fmin = minimize(likelihood1, line, model, wave, flux, ivar, noise, x, **init_pars2)

    return m, fmin

def plot_fit(wave,line, model, noise, m, color, x=None, lab=None):
    plt.plot(wave, line(*[m.values[p] for p in m.parameters], wave=wave, x=x, model=model, noise= noise), color, lw=2, alpha = .8, label = lab)
    
def get_chi(like, line, model, m, wave, flux, ivar, noise = None, x=None):
    return str(like(*[m.values[p] for p in m.parameters], line = line, model=model, noise= noise, wave=wave, flux = flux, ivar = ivar, x = x)) +' / ' + str((len(wave) - len(m.parameters)))

def rebin(x, wind, wave, flux):
    A = wave - x[:,None]
    w = (A >= -wind/2) & (A < wind/2) # np.abs(A < dlam/2) #
    A[w] = 1
    A[~w] = 0
    norma  = A.sum(axis = 1) # norm for each row
    norm_nul = (norma == 0) #create mask where the nrom is 0
    norma[norm_nul] = 1 # if norm ==0 set it to 1 for division
    A = A / norma[:,None]
    flux = A.dot(flux)
    return flux

def init_model(wave, flux, window, model, model_label):
    par_names = model.parnames
    if model_label == 'spl':
        pos_max = wave[np.where(flux == flux.max())[0][0]]
        x_node = np.arange(pos_max - window, pos_max + window, 2*window/len(par_names))
        if len(x_node) != len(par_names):
            x_node = x_node[1:]
        tik = rebin(x_node, 10, wave, flux)
        init_pars = {x:y for x,y in zip(par_names, tik)}
        plt.plot(x_node, tik,'.k')
    else:
        init_pars = {x:y for x,y in zip(par_names, [1, wave.mean(), 10, 0, 0])}
        x_node = None
        if model_label == 'asym_lorentzian':
            init_pars['c2'] = 10
    return x_node, init_pars
    
