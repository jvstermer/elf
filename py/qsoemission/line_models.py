import numpy as np
    
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
    # a : amplitude of the gaussian
    # b : mean of the lorentzian
    # c : width of the lorentzian
    # d : slope of the linear function
    # e : intercept of the linear function
    # returns lorentzian + linear 
    
    wave = kwargs['wave']
    return a / np.pi * .5 * c / ( (wave - b)**2 + (c/2)**2 ) + d * wave + e

