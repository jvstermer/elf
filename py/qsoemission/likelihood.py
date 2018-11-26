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
    
    line_model = kwargs['line_model']
    wave = kwargs['wave']
    flux = kwargs['flux']
    ivar = kwargs['ivar']
    
    res = (flux-line_model(*args, wave=wave))*np.sqrt( ivar )
    return res.dot(res.T)
    
def like_new(*args, **kwargs):
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

    line_model = kwargs['line_model']
    ivar = kwargs['ivar']
    
    mask = ivar != 0
    ivar = ivar[mask]
    
    wave = kwargs['wave'][mask]
    flux = kwargs['flux'][mask]
    
    v = line_model(*args, wave = wave)

    y = sps.erfc(flux / ( np.sqrt(2) * ivar)) - sps.erfc((flux - v) / ( np.sqrt(2) * ivar))
    
    mask1 = y != 0
    
    return np.sum(np.log(ivar[mask1]) + np.log(np.abs(v[mask1])) - np.log(np.abs(y[mask1])))
    #return np.sum( np.log(np.abs(v[mask1])) - np.log(np.abs(y[mask1])))

