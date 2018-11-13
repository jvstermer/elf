import numpy as np

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
    
def vrai(*args, **kwargs):

    line_model = kwargs['line_model']
    wave = kwargs['wave']
    flux = kwargs['flux']
    ivar = kwargs['ivar']
    
    v = line_model(*args, wave = wave)
    
    x1 = flux / ( np.sqrt(2) * ivar)
    x2 = (flux - v) / ( np.sqrt(2) * ivar)
    y = sps.erfc(x1) - sps.erfc(x2)
    
    dist = np.sqrt( np.pi / 2 ) * ivar * y / v
    return -np.log(np.prod(dist))

