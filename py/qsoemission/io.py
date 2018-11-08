import numpy as np
import glob
import fitsio

def read_pix(filename, qso_id):
    h = fitsio.FITS(filename)
    flux =  h['FLUX'][:,qso_id].T[0]
    wave = 10**h['LOGLAM_MAP'][:] # in LF
    ivar = h['IVAR'][:,qso_id].T[0]
    qso_num = h['THING_ID_MAP'][:][qso_id]
    return flux, ivar, wave, qso_num

def read_drq(drq):
    vac = fitsio.FITS(drq)
    try:
        z_dict = {x:y for x,y in zip(vac[1]['THING_ID'][:], vac[1]['Z'][:])}
    except:
        z_dict = {x:y for x,y in zip(vac[1]['THING_ID'][:], vac[1]['Z_VI'][:])}    
    return z_dict
    
