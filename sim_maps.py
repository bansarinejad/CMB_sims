import numpy



def create_correlated_power(cov_4d):
    '''Will be used by both the noise and FG codes
    
    dimensions of cov_4d are npad x npad x nf x nf
    
    '''

    npad = cov_4d.shape[0]
    nfreq = cov_4d.shape[2]
    gaussian_real = np.random.normal(size = [npad,npad,nfreq])

        #now a bunch of matrix multiplies
    Cholesky = np.linalg.cholesky(cov_4d)
    out_maps = np.einsum(Cholesky, gaussian_real, stuff to define indices)
    return out_maps


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
        self.nfreq=nfreq
        self.ncombo = nfreq * (nfreq+1)/2
        self.pairs = np.zeros([self.ncombo,2],dtype=int)
        k=0
        for i in range(nfreq):
            for j in range(i,nfreq):
                self.pairs[k,0]=i
                self.pairs[k,1]=j
                k+=1
        
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
    
    def create_cmb(self,Cov_4d = None):
        '''
        what basis? if IQU will need to change this code to iterate over 3-vectors, which is doing T/E (2-vector) here
        '''
        if Cov_4d is None:
            Cov_4d = np.zeros([self.npad,self.npad,2,2])
            theory_psd = self.camb_to_psd(0,0)
            Cov_4d[:,:,0,0]=theory_psd
            theory_psd = self.camb_to_psd(1,1)
            Cov_4d[:,:,1,1]=theory_psd
            theory_psd = self.camb_to_psd(0,1)
            Cov_4d[:,:,0,1]=Cov_4d[:,:,1,0]=theory_psd

        return create_correlated_power(Cov_4d)

    def create_fgs(self,Cov_4d = None):
        '''
        Assumes have a set of Theory Spectra
        Will turrn each spectrum to a PSD that is npad x npad in size
        and combine these to create correlated noise
        eg 90x90;, 90x150, 90x220, 150x150, 150x220, 220x220)
        ordering from self.pairs

        Can also pass in an existing Cov structure, in which case, it skips to creating the realization. 
        '''
        if Cov_4d is None:
            Cov_4d = np.zeros([self.npad,self.npad,self.nfreq,self.nfreq])
            k=0
            for i in range(self.nfreq):
                for j in range(i,self.nfreq):
                    theory_psd = self.turn_spectra_to_psd(i,j)
                    Cov_4d[:,:,i,j] = Cov_4d[:,:,j,i] = theoy_psd
                    k+=1

        return create_correlated_power(Cov_4d)


    def create_noise(self,Cov_4d = None):
        '''
        Assumes have a set of PSDs
        An individual PSD is npad x npad in size
        For 3 observing freqs, have 6 PSDs --
        eg 90x90;, 90x150, 90x220, 150x150, 150x220, 220x220)
        ordering from self.pairs
        '''
        if Cov_4d is None:
            Cov_4d = np.zeros([self.npad,self.npad,self.nfreq,self.nfreq])
            k=0
            for i in range(self.nfreq):
                for j in range(i,self.nfreq):
                    Cov_4d[:,:,i,j] = Cov_4d[:,:,j,i] = self.psds[k,:,:]
                    k+=1

        return create_correlated_power(Cov_4d)

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
