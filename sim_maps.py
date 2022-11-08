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

        reso = self.reso_arcmin * np.pi / 180 / 60.
        self.dell_pad = 2*np.pi/(self.npad * reso)
        self.ell1d = 2*np.pi * np.fft.fftfreq(self.npad,d=reso)
    
    def create_cmb(self):
        pass

    def create_fgs(self):
        pass

    def create_noise(self):
        pass

    def lens_maps(self, large_maps):
        pass

    def filter_map(self, large_map, index = None, map_in_fourier = True):
        '''Apply filter in Fourier space. Does not return to real space'''
        if not map_in_fourier:
            raise Exception("Haven't implemented real map in filter_maps")
        assert large_map.shape[0] == self.npad
        assert large_map.shape[1] == self.npad
        if index is None:
            return large_map * self.filters
        return large_map * self.filters[index,:,:]
        

    def cut_map(self, large_map):
        '''Cut maps in real space from npad x npad, taking one corner of nfinal x nfinal'''
        return large_map[:self.nfinal,:self.nfinal]
 

    def zeropad_map(self,map):
        ''' Zero pad map'''
        large_map =  np.zeros([self.npad,self.npad])

        large_map[:self.nfinal,:self.nfinal] = map
        pass

    def apodize_map(self,map):
        '''Mean subtract and apodize small map'''
        mn = np.average(map)
        out_map = np.zeros([self.nfinal,self.nfinal])
        out_map = map - mn
        out_map *= self.apod
        return out_map

    def zeropad_and_apodize_map(self,map):
        ''' Meansubtract, apodize and zeropad map'''
        return zeropad_map(apodize_map(map))

    def source_interpolate(self, map, source_mask):
        pass

    def apply_beams(self,large_maps, index = None, map_in_fourier = True):
        '''Apply beams in Fourier space. Does not return to real space'''
        if not map_in_fourier:
            raise Exception("Haven't implemented real map in apply_beams")
        assert large.maps.shape[0] == self.npad
        assert large.maps.shape[1] == self.npad
        lx = np.tile(self.ell1d, [1,self.npad])
        #should be npad x npad now
        ll = np.sqrt(lx*lx  + lx.T*lx.T)
        lbl = self.lbeam
        if beam_index is None:
            bl = self.beams[:]
        else:
            bl = self.beams[index,:]
        bl2d = np.interp(ll, lbl,bl,left=0,right=0)
        return large_maps * bl2d
