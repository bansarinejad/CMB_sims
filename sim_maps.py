import numpy


class patch:
    def __init__(self,reso_arcmin=0.5, npad = 256, nfinal = 128, nfreq = 3, lmax = 15000,
                camb_file = None, fg_file = None, psd_file = None, filter_file = None,apod_file=None, beam_file=None) -> None:
        self.reso_arcmin = reso_arcmin
        self.npad = npad
        self.nfinal = nfinal
        self.camb_file = camb_file
        self.fg_file = fg_file
        self.psd_file = psd_file
        self.filter_file = filter_file
        self.apod_file = apod_file
        self.beam_file = beam_file
        self.lmax = lmax
        #splaceholders for loading spectra
        self.fg_spectra = load_fgs(fg_file,lmax)
        self.cmb_spectra = load_camb(camb_file, lmax)
        #palceholders for loading PSD and filters:
        #npad/reso passed in for error checking against what was used in creating the files
        self.psds = load_psds(psd_file, lmax, npad, reso_arcmin)
        self.filters = load_filters(psd_file, lmax, npad, reso_arcmin)
        self.apod = load_apod(apod_file)
        self.beams = load_beams(beam_file)
    
    def create_cmb(self):
        pass

    def create_fgs(self):
        pass

    def create_noise(self):
        pass

    def lens_maps(self, large_maps):
        pass

    def filter_maps(self, large_maps):
        pass

    def cut_maps(self, large_maps):
        pass

    def zeropad_maps(self,maps):
        pass

    def apodize_maps(self,large_maps):
        pass

    def source_interpolate(self, maps, source_mask):
        pass

    def apply_beams(self,large_maps):
        pass