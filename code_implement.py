import numpy as np
import matplotlib.pyplot as plt
import const
import iminuit
import read_file
import line_models
import likelihood
from configparser import ConfigParser

loc = "./util.ini"

def go_to_lf(lam, z):
    return np.log10( lam * (z+1) )

def set_system_values(path, section, value):

    cp = ConfigParser()       
    cp.read(path)
    return cp.get(section,value)
    

def window(z, flux, wave, ivar, line_id):

    range_window = float(set_system_values(loc, 'chi2 scan', 'window'))
    line = const.emission_lines[line_id]
    mask = ( wave >= go_to_lf(line - range_window /2, z) ) & ( wave <= go_to_lf(line + range_window /2, z) )
    mi = minimization(wave[mask], flux[mask], ivar[mask], go_to_lf(line, z))
    plt.plot(wave[mask], mi)
    return mi, 10** mi[1] / line - 1
    

def minimization(x, flux, ivar, line_lam):

    def model(a,b,c,d,e):
        return getattr(line_models, set_system_values(loc,'line model','model'))
    
    m = iminuit.Minuit(
        getattr(likelihood, set_system_values(loc, 'likelihood method', 'method')),
        b = line_lam,
        c = float(set_system_values(loc, 'minimization', 'sigma')),
        pedantic = False,
        )
        
    fmin  = m.migrad()
    
    inte = list(m.values.values()) 

    return model(inte)

class qso:

    def __init__(self, file_id, qso_id):
        pix_file_list = set_system_values(loc, 'input','filename')
        self.flux, self.ivar, self.wave, self.id = read_file.read_pix(pix_file_list, file_id, qso_id)
        return

file_id = 0
qso_id = 480
qso1= qso(file_id, qso_id)

z_dict = read_file.read_drq(set_system_values(loc, 'input','z_file'))
z = z_dict[qso1.id]

if qso1.id in z_dict and z > 2.1:
    plt.plot(qso1.wave, qso1.flux)
    fit, new_z = window( z, qso1.flux, qso1.wave, qso1.ivar,'Lya')
    plt.show()

# define wavelength range where to fit the emission line
    # z : global redshift from file
    # flux : qso flux
    # wave : qso wavelength range
    # line_id : which emission line to fit
    # returns : result of the fit and the new redshift
