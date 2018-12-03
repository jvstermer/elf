import numpy as np
import scipy.special as sps


def chi_squared(*args, **kwargs):
    '''
    Calculates the chi2 between the flux and the model
    using ivar as inverse variances

    Arguments:
        wave -- wavelength (numpy array)
        flux -- flux (numpy array)
        ivar -- inverse variance (numpy array)
        line_model -- the model for the line (function)
        **kwargs -- keyword arguments to be passed to line_model

    Returns:
        the chi2 (float)
    '''    
    line = kwargs['line']
    flux = kwargs['flux']
    ivar = kwargs['ivar']
    
    res = (flux-line(*args, **kwargs))*np.sqrt( ivar )
        
    return res.dot(res.T)
    
def new_like(*args, **kwargs):
    '''
    Calculates the likelihood of the model
    using ivar as inverse variances

    Arguments:
        wave -- wavelength (numpy array)
        flux -- flux (numpy array)
        ivar -- inverse variance (numpy array)
        line_model -- the model for the line (function)
        **kwargs -- keyword arguments to be passed to line_model

    Returns:
        the likelihood (float)
    '''
    
        
    line = kwargs['line']
    ivar = kwargs['ivar']
    
    mask = ivar != 0
    ivar = ivar[mask]
    flux = kwargs['flux'][mask]

    v = line(*args,**kwargs)
    v = v[mask]
    y = sps.erfc(flux / ( np.sqrt(2) * ivar)) - sps.erfc((flux - v) / ( np.sqrt(2) * ivar))
    
    mask1 = y != 0
    
    return np.sum(np.log(ivar[mask1]) + np.log(np.abs(v[mask1])) - np.log(np.abs(y[mask1])))
    

