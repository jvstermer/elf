emission_lines = {
    'Lya'   : 1215.24,
    'CIII'  : 1908.734,
    'CIV'   : 1549.48,
    'MgII'  : 2799.117,
    }

line_dict = {
    'Lya_z' : [],
    'Lya_err' : [],
    'CIII_z' : [],
    'CIII_err' : [],
    'CIV_z' : [],
    'CIV_err' : [],
    'MgII_z' : [],
    'MgII_err' : [],
    'cat_z' : [],
}

c = 299792.458 # km/s

# number of parameters for set functions
dict_num_pars = {
    'gaussian'          : 3, 
    'lorentzian'        : 3, 
    'asym_lorentzian'   : 4, 
    }
