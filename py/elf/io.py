from __future__ import print_function
import numpy as np
import scipy as sp
import glob
import fitsio
import healpy
import sys
import time
import os.path

from elf import const, qso, utils

###################################################
######################  READ  #####################
###################################################

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

def read_new_spplate(h, which_qso):
    ivar = h[1]['ivar'][which_qso]
    flux = h[1]['flux'][which_qso]
    wave = h[1]['loglam'][which_qso]
    z = h[1]['z'][which_qso]
    return flux, ivar, wave, z

def dict_drq(drq):
    vac = fitsio.FITS(drq)
    try:
        z_dict = {x:y for x,y in zip(vac[1]['THING_ID'][:], vac[1]['Z'][:])}
    except:
        z_dict = {x:y for x,y in zip(vac[1]['THING_ID'][:], vac[1]['Z_VI'][:])}    
    return z_dict

def read_drq(drq, zmin):
    vac = fitsio.FITS(drq)

    ## Redshift
    try:
        zqso = vac[1]["Z"][:]
    except:
        sys.stderr.write("Z not found (new DRQ >= DRQ14 style), using Z_VI (DRQ <= DRQ12)\n")
        zqso = vac[1]["Z_VI"][:]

    ## Info of the primary observation
    thid  = vac[1]["THING_ID"][:]
    ra    = vac[1]["RA"][:].astype('float64')
    dec   = vac[1]["DEC"][:].astype('float64')
    plate = vac[1]["PLATE"][:]
    mjd   = vac[1]["MJD"][:]
    fid   = vac[1]["FIBERID"][:]
    
    w = (thid>0)
    w &= (ra!=dec)
    w &= (ra!=0.)
    w &= (dec!=0.)
    w &= (zqso>0.)

    ## Redshift range
    if not zmin is None:
        w &= zqso>zmin

    ra    = ra[w]*sp.pi/180.
    dec   = dec[w]*sp.pi/180.
    zqso  = zqso[w]
    thid  = thid[w]
    plate = plate[w]
    mjd   = mjd[w]
    fid   = fid[w]
    vac.close()

    return ra,dec,zqso,thid,plate,mjd,fid

def read_data(in_dir, drq, zmin = None, log = False, start_plate = None, end_plate = None):

    ra,dec,zqso,thid,plate,mjd,fid = read_drq(drq,zmin)

    data = read_from_spplate(in_dir,thid, ra, dec, zqso, plate, mjd, fid, start_plate, end_plate)

    flux = [d.fl for d in data] ; flux = sp.array(flux)
    ivar = [d.iv for d in data] ; ivar = sp.array(ivar)
    wave = [d.ll for d in data] ; wave = sp.array(wave)
    
    if log == False:
        wave = 10**wave
    
    if data == []:
        zqso = []
        
    return flux, ivar, wave, zqso

def read_from_spplate(in_dir, thid, ra, dec, zqso, plate, mjd, fid, start_plate, end_plate):
    pix_data={}
    unique_plates = sp.unique(plate)[start_plate:end_plate]

    for i,p in enumerate(unique_plates):
        sys.stdout.write('\rplate %d: %d / %d '%(p, i+1, len(unique_plates)))
        sys.stdout.flush()
        
        wplate = plate==p

        spplates = glob.glob(in_dir+"/{}/spPlate-{}-*.fits".format(p, p))

        for spplate in spplates:
            h = fitsio.FITS(spplate)
            head0 = h[0].read_header()
            MJD = head0["MJD"]
            
            wfib = wplate
            wmjd = mjd == MJD
            wfib = wplate & wmjd
            
            coeff0 = head0["COEFF0"]
            coeff1 = head0["COEFF1"]

            flux = h[0].read()
            ivar = h[1].read()*(h[2].read()==0) #sets all masked ivar to 0
            llam = coeff0 + coeff1*sp.arange(flux.shape[1])
            
            for (t, r, d, z, p, m, f) in zip(thid[wfib], ra[wfib], dec[wfib], zqso[wfib], plate[wfib], mjd[wfib], fid[wfib]):
                index = f-1
                d = qso.forest(llam,flux[index],ivar[index])
                if t in pix_data:
                    pix_data[t] += d
                else:
                    pix_data[t] = d
            h.close()

    return list(pix_data.values())

###################################################
######################  WRITE  ####################
###################################################

def write_in_dict(model, dici, l, m, wa, x_node= None ):
    if model.label == 'spl':
        y = model(*[m.values[p] for p in m.parameters], wave=wa, x=x_node)
        amp = y.max()
        dici[l+'_z'].append(wa[np.where(y == amp)[0][0]] / const.emission_lines[l] - 1)
        dici[l+'_err'].append(-1)#amp / m.errors['a'])     
    else:   
        dici[l+'_z'].append(m.values['a1'] / const.emission_lines[l] - 1)
        if m.errors['a0'] == 0:
            m.errors['a0'] = 1
        dici[l+'_err'].append(m.values['a0'] / m.errors['a0'])
        
    return dici

def write_fits(out, dic):
    fits = fitsio.FITS(out,'rw', clobber = True)
    
    na = list(dic.keys())
    lis = [np.array(dic[key])  for key in dic.keys()]
    fits.write(lis, names = na)
    fits.close()

