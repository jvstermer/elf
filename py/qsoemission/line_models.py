import numpy as np

def linear( d, e, **kwargs):

    wave = kwargs['wave']
    return d * wave + e
    
def quad(d,e,f, **kwargs):

    wave = kwargs['wave']
    return f * wave**2 + d * wave + e
    
    
def gaussian_lin(a, b, c, d, e, **kwargs):
    #takes array of values
    # a : amplitude of the gaussian
    # b : mean of the gaussian
    # c : standard deviation of the gaussian
    # d : slope of the linear function
    # e : intercept of the linear function
    # returns gaussaian + linear 

    wave = kwargs['wave']
    return a*np.exp(-(wave-b)**2/2/c**2 ) + d * wave + e

def lorentzian_lin(a, b, c, d, e, **kwargs):
    #takes array of values
    # a : amplitude of the lorentzian
    # b : mean of the lorentzian
    # c : width of the lorentzian
    # d : slope of the linear function
    # e : intercept of the linear function
    # returns lorentzian + linear 
    
    wave = kwargs['wave']
    #return a  * c / 2 / ( (wave - b)**2 + (c/2)**2 ) + d * wave + e
    return a / (1+((wave-b)*2 /c)**2) + d * wave + e

def asym_lorentzian(a, b, c, c2, d, e, **kwargs):
    #takes array of values
    # a : amplitude of the lorentzian
    # b : mean of the lorentzian
    # c1 : lefthand width of the lorentzian
    # c2 : rightthand width of the lorentzian
    # d : slope of the linear function
    # e : intercept of the linear function
    # returns asymmetrical lorentzian + linear 
    
    wave = kwargs['wave']
    
    f = lorentzian_lin(a, b, c, d, e, wave = wave)
    f[np.where(wave>b)] =  lorentzian_lin(a, b, c2, d, e, wave = wave[np.where(wave>b)])
    return f

def lorentzian_quad(a, b, c, d, e, f, **kwargs):
    #takes array of values
    # a : amplitude of the lorentzian
    # b : mean of the lorentzian
    # c : width of the lorentzian
    # d : slope of the linear function
    # e : intercept of the linear function
    # returns lorentzian + linear 
    
    wave = kwargs['wave']
    return a / np.pi * .5 *c / ( (wave - b)**2 + (c/2)**2 ) + f * wave**2 + d * wave + e
    
    
    
