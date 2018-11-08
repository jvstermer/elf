import numpy as np

def chi_squared(wave=None, flux=None, ivar=None, line_model=None, *args):
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
    
    res = (flux-line_model(wave=wave, **kwargs))*np.sqrt( ivar )
    return res.dot(res.T)
# add new method
