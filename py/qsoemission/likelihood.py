import numpy as np

"""def chi_squared(flux, ivar, model, a,b,c,d,e):
    res = ( flux - model (a,b,c,d,e) ) * np.sqrt( ivar )
    return res.dot(res.T)"""
def chi_squared(flux, ivar, model, co):
    res = ( flux - model (co) ) * np.sqrt( ivar )
    return res.dot(res.T)
# add new method
