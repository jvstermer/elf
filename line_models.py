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
    
def gaussian_lin(x,co):
    #takes array of values
    # a : amplitude of the gaussian
    # b : mean of the gaussian
    # c : standard deviation of the gaussian
    # d : slope of the linear function
    # e : intercept of the linear function
    # returns gaussaian + linear 
    return co[0] * np.exp(-(x-co[1])**2 / 2 * co[2]**2 ) + co[3] * x + co[4]
