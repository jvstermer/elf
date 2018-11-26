import numpy as np
import glob
import fitsio
from py.qsoemission import const

def read_delta(j):
    ivar = j['WEIGHT'][:]
    flux =  j['DELTA'][:]
    wave = 10**j['LOGLAM'][:] # in LF
    head = j.read_header()
    z = head["Z"]
    return ivar, flux, wave, z

def read_pix(h, qso_id):
    ivar = h['IVAR'][:,qso_id].T[0]
    flux =  h['FLUX'][:,qso_id].T[0]
    wave = 10**h['LOGLAM_MAP'][:] # in LF
    qso_num = h['THING_ID_MAP'][:][qso_id]
    return flux, ivar, wave, qso_num

def read_drq(drq):
    vac = fitsio.FITS(drq)
    try:
        z_dict = {x:y for x,y in zip(vac[1]['THING_ID'][:], vac[1]['Z'][:])}
    except:
        z_dict = {x:y for x,y in zip(vac[1]['THING_ID'][:], vac[1]['Z_VI'][:])}    
    return z_dict

def write_in_dict(dici, l, m ):
    dici[l+'_z'].append(m.values['b'] / const.emission_lines[l] - 1)
    
    if m.errors['a'] == 0:
        m.errors['a'] = 1
    dici[l+'_err'].append(m.values['a'] / m.errors['a'])
    
    """if m.errors['a'] == 0:
        dici[l+'_err'].append(m.values['a'])
    else:
        dici[l+'_err'].append(m.values['a'] / m.errors['a'])"""
    return dici

def write_fits(out, dic):
    fits = fitsio.FITS(out,'rw', clobber = True)
    
    na = list(dic.keys())
    lis = [np.array(dic[key])  for key in dic.keys()]
    fits.write(lis, names = na)
    fits.close()
