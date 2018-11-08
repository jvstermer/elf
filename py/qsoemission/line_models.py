import numpy as np

"""def gaussian_lin(x,a,b,c,d,e):
    #takes array of values
    # a : amplitude of the gaussian
    # b : mean of the gaussian
    # c : standard deviation of the gaussian
    # d : slope of the linear function
    # e : intercept of the linear function
    # returns gaussaian + linear 
    return a * np.exp(-(x-b)**2 / 2 * c**2 ) + d * x + e"""
    
def gaussian_lin(wave, **kwargs):
    #takes array of values
    # a : amplitude of the gaussian
    # b : mean of the gaussian
    # c : standard deviation of the gaussian
    # d : slope of the linear function
    # e : intercept of the linear function
    # returns gaussaian + linear 

    a = kwargs['a']
    b = kwargs['b']
    c = kwargs['c']
    d = kwargs['d']
    e = kwargs['e']

    return a*np.exp(-(wave-b)**2/2/c**2 ) + d*wave + e
