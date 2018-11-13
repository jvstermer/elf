#!/usr/bin/env python
from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sps
import iminuit
from configparser import ConfigParser
import sys
import glob
import lib
from py.qsoemission import io, utils, qso, line_models, likelihood, const

if len(sys.argv) != 2:
    print('usage: find_lines config.ini')
    sys.exit(0)

z_file = utils.get_system_values(sys.argv[1], 'input','z_file')
z_dict = io.read_drq(z_file)

window = utils.get_system_values(sys.argv[1], 'chi2 scan', 'window')
window = float(window)

filename = glob.glob(utils.get_system_values(sys.argv[1], 'input', 'filename'))[int(utils.get_system_values(sys.argv[1], 'input', 'pix_file_id'))]
qso_id = int(utils.get_system_values(sys.argv[1], 'input','qso_id'))
qso1= qso.qso(filename, qso_id)
z = z_dict[qso1.id]

def line_model(a, b, c, d, e, wave):
    return a * np.exp(-(wave-b)**2 / 2 * c**2 ) + d * wave + e

def vrai(a, b, c, d, e, line_model, wave, flux, ivar):
    
    v = line_model(a, b, c, d, e, wave)
    
    x1 = flux / ( np.sqrt(2) * ivar)
    x2 = (flux - line_model(a, b, c, d, e, wave)) / ( np.sqrt(2) * ivar)
    y = x1 - x2
    
    dist = np.sqrt( np.pi / 2 ) * ivar * y / v
    return -np.log(np.prod(dist))


l = 'Lya'
wa, fl, iv = utils.window(z, qso1.wave, qso1.flux, qso1.ivar, l, window)
#a = 25.785997251685714 ; b = 4307.586484952259 ; c = 12.5319602614264 ; d = -0.0001553733277506994 ; e = 4.868459789227105 
#print(vrai(a, b, c, d, e, line_model, wa, fl, iv))
m, fmin = utils.minimize(vrai, line_model, wa, fl, iv, a=25, b=wa.mean(), c=10., d=0, e=5, pedantic = False)

