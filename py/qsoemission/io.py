import numpy as np
import glob
import fitsio

def read_pix(init, file_id, qso_id):
    to_open = glob.glob(init)
    h = fitsio.FITS(to_open[file_id])
    flux =  h['FLUX'][:,qso_id].T[0]
    wave = h['LOGLAM_MAP'][:] # in LF
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
    
