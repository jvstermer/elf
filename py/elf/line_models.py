import numpy as np
import scipy.interpolate as interpolate

class line_model:
    def __init__(self, func, parnames, lab):
        self.parnames = parnames
        self.func = func
        self.label = lab
        
    def __call__(self, *pars, **kwargs):
        return self.func(*pars, **kwargs)

def add(*args, **kwargs ):
    model = kwargs['model']
    noise = kwargs['noise']
    
    x = model(*[p for p in args[:len(model.parnames)]], **kwargs)
    if noise != None:
        x += noise(*[p for p in args[len(model.parnames):]], **kwargs)
    return x

def polynomial(*args, **kwargs):
    wave = kwargs['wave']
    
    res = 0
    for index, coeff in enumerate(args):
        res += coeff * wave** index
    return res 

def gaussian(*pars,**kwargs):
    #takes array of values
    # a : amplitude of the gaussian
    # b : mean of the gaussian
    # c : standard deviation of the gaussian
    # returns gaussaian

    wave = kwargs['wave']
    return pars[0]*np.exp(-(wave-pars[1])**2/2/pars[2]**2 )

def lorentzian(*args, **kwargs):
    #takes array of values
    # a : amplitude of the lorentzian
    # b : mean of the lorentzian
    # c : width of the lorentzian
    # returns lorentzian
    
    wave = kwargs['wave']
    return args[0] / (1+((wave-args[1])*2 /args[2])**2)

def asym_lorentzian(*args, **kwargs):
    #takes array of values
    # a : amplitude of the lorentzian
    # b : mean of the lorentzian
    # c1 : lefthand width of the lorentzian
    # c2 : rightthand width of the lorentzian
    # returns asymmetrical lorentzian 
    
    wave = kwargs['wave']
    
    f = lorentzian(args[0], args[1], args[2],wave = wave)
    f[np.where(wave>args[1])] =  lorentzian(args[0], args[1], args[3], wave = wave[np.where(wave>args[1])])
    return f
        
def spl(*pars, **kwargs):
        wave = kwargs['wave']
        x = kwargs['x']
        return interpolate.splev( wave, interpolate.splrep(x, pars, s=0, k=3))
        

