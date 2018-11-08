import matplotlib.pyplot as plt ; import numpy as np ; import scipy as sp
import fitsio ; import glob ; import iminuit

"""def minimization(x,flux, line, sig, ivar):"""

    #def dist(a):
        #return 1

    def vrai(c, sig):
        m = np.zeros(len(flux))
        def integrand(flux, c, a, sig):
            return np.exp( - ( flux - c * a ) **2 / ( 2 * sig **2 ) ) # * dist(a)
        for i in range(len(flux)):
            m[i], _ = quad( integrand, 0, 1, args = (flux[i], c[i], sig) )
        return np.prod(m)

    """mig = iminuit.Minuit(vrai, pedantic = False)
    mig.print_param()
    fmin = mig.migrad()
    inte = list(mig.values.values())
    print(inte)
    return vrai(inte),inte[1]"""
    
#################Where to do the fit and plot it######################### 
def window(z, flux, wave, line, sigma, range_window,ivar):
    mask = ( wave >= np.log10((line-range_window/2) * (z+1)) ) & ( wave <=np.log10((line+range_window/2) * (z+1)) ) #RF -> LF
    
    mi, b = minimization(wave[mask], flux[mask], np.log10(line * (z+1)), sigma, ivar[mask])
    new = 10**b / line - 1
    plt.plot(wave[mask], mi)

    return new

##############get data#####################
list_of_filenames = glob.glob( '../original_files/pix/pix_*.fits' )
drq = fitsio.FITS('../original_files/DR14Q_v4_4.fits')

######DRQ############
z_dict = {x:y for x,y in zip(drq[1]['THING_ID'][:], drq[1]['Z'][:])}

names = ['La','MgII', 'CIII', 'CIV' ]
lines = [1215.24, 2799.117, 1908.734, 1549.48]

h = fitsio.FITS(list_of_filenames[0])
j = 480
z = z_dict[h['THING_ID_MAP'][:][j]]

new_z = []

if h['THING_ID_MAP'][:][j] in z_dict and z > 2.1:
    flux =  h['FLUX'][:,j].T[0]
    wave = h['LOGLAM_MAP'][:] # in LF
    ivar = h['IVAR'][:,j].T[0]
    plt.plot(wave,flux, alpha = .5)
    red = window(z, flux, wave, lines[0], 3, 200, ivar)
    """for i in range(len(lines)):
        red = window(z, flux, wave, lines[i], 3, 200, ivar)
        new_z.append(red)
        print (chi/(card - 5))"""
    ##########Plot############################################
    #plt.show()

