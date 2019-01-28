from . import io

class qso:
    def __init__(self, filename, which_qso, mode):
        if mode == 'pix':
            self.flux, self.ivar, self.wave, self.id = io.read_pix(filename, which_qso)
        if mode == 'spplate':
            self.flux, self.ivar, self.wave, self.z = io.read_new_spplate(filename, which_qso)

class forest:
    def __init__(self,ll,fl,iv):
        self.ll = ll
        self.fl = fl
        self.iv = iv

