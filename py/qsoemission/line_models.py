import numpy as np
import scipy.interpolate as interpolate

class line_model:
    def __init__(self, func, parnames):
        self.parnames = parnames
        self.func = func
        
    def __call__(self, *pars, **kwargs):
        return self.func(*pars, **kwargs)

def add(*args, **kwargs ):
    model = kwargs['model']
    noise = kwargs['noise']
    
    x = model(*[p for p in args[:len(model.parnames)]], **kwargs)
    if noise != None:
        x += noise(*[p for p in args[len(model.parnames):]], **kwargs)
    return x

def linear( d, e, **kwargs):
    wave = kwargs['wave']
    return d * wave + e

def gaussian(a, b, c,**kwargs):
    #takes array of values
    # a : amplitude of the gaussian
    # b : mean of the gaussian
    # c : standard deviation of the gaussian
    # returns gaussaian

    wave = kwargs['wave']
    return a*np.exp(-(wave-b)**2/2/c**2 )

def lorentzian(a, b, c, **kwargs):
    #takes array of values
    # a : amplitude of the lorentzian
    # b : mean of the lorentzian
    # c : width of the lorentzian
    # returns lorentzian
    
    wave = kwargs['wave']
    return a / (1+((wave-b)*2 /c)**2)

def asym_lorentzian(a, b, c, c2, **kwargs):
    #takes array of values
    # a : amplitude of the lorentzian
    # b : mean of the lorentzian
    # c1 : lefthand width of the lorentzian
    # c2 : rightthand width of the lorentzian
    # returns asymmetrical lorentzian 
    
    wave = kwargs['wave']
    
    f = lorentzian(a, b, c,wave = wave)
    f[np.where(wave>b)] =  lorentzian(a, b, c2, wave = wave[np.where(wave>b)])
    return f
        
def spl(*pars, **kwargs):
        wave = kwargs['wave']
        x = kwargs['x']
        return interpolate.splev( wave, interpolate.splrep(x, pars, s=0, k=3))
