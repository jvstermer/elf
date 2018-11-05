import numpy as np

def chi_squared(a,b,c,d,e):
    res = ( flux - model (a,b,c,d,e) ) * np.sqrt( ivar )
    return res.dot(res.T)
# add new method
