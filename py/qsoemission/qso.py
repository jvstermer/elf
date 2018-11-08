from . import io

class qso:
    def __init__(self, filename, qso_id):
        self.flux, self.ivar, self.wave, self.id = io.read_pix(filename, qso_id)
