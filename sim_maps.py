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
    #Cholesky is i, j, k, l
    #gauss_real is i, j, l
    out_maps = np.einsum('ijkl,ijl->ijk',Cholesky, gaussian_real)
    return out_maps

def fill_in_beams_fromfile(files):
    nl = len(ells)
    nfreq = len(files)
    beams_interp = np.zeros([nspectra, nl],dtype=np.float32)
    for i in range(nfreq):
        bl = np.fromfile(files[i]) # this will need to be fixed
        beams_interp[i,:] = np.interp(ells,bl[0,:],bl[1,:])
        bad = beams_interp[i,:] < 0
    return beams_interp

def fill_in_theory(files,ells,cl2dl=False):
    nl = len(ells)
    
    nfreq = len(files)
    theory_interp = np.zeros([nfreq, nl],dtype=np.float32)
    for i in range(nfreq):
        dls = np.loadtxt(files[i])
        locl=dls[:,0]
        dl = dls[:,1]
        if cl2dl:
            dl = dl * locl*(locl+1)/(2*np.pi)
        dl[locl < 2]=0 #don't want 0,1
        theory_interp[i,:] = np.interp(ells,locl,dl)
    return theory_interp   

class patch:
    def __init__(self,reso_arcmin=0.5, npad = 256, nfinal = 128, nfreq = 3, lmax = 15000,
                 freq_names=['95,'150','220'],
                camb_file = None, fg_file = None, psd_file = None, filter_file = None,apod_file=None, 
                beam_file=None) -> None:
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
        if lmax >= 15000:
            print('Warning: default spt3g beam file stopped at l=14999, vs requested lmax = {}'.format(lmax))
        self.ells = range(lmax+1)
        self.nfreq=nfreq
        self.freq_names=freq_names
        self.ncombo = nfreq * (nfreq+1)/2
        self.pairs = np.zeros([self.ncombo,2],dtype=int)
        k=0
        for i in range(nfreq):
            for j in range(i,nfreq):
                self.pairs[k,0]=i
                self.pairs[k,1]=j
                k+=1
        
        #splaceholders for loading spectra
        self.fg_spectra = (np.loadtxt(fg_file)).T #this is expected to yield an array: [1+nfreq, Nl]. 0 is l, 1-Ncombo is 90x90, 90x150, etc.
        self.cmb_spectra = (np.loadtxt(camb_file)).T #this is expected to yield an array: [1+nfreq, Nl]. 0 is l, 1 is TT, 2 is EE, 3 is BB, 4 is TE
        #palceholders for loading PSD and filters:
        #npad/reso passed in for error checking against what was used in creating the files
        self.psds = load_psds(psd_file, lmax, npad, reso_arcmin)
        self.filters = load_filters(psd_file, lmax, npad, reso_arcmin)
        self.apod = load_apod(apod_file)
        # expect something like: '/home/creichardt/spt3g_software/beams/products/compiled_2020_beams.txt'
        self.beams = (np.loadtxt(beam_file)).T  #this is expected to yield an array: [1+nfreq, Nl]. Note the above file's last ell = 14999
        

        reso = self.reso_arcmin * np.pi / 180 / 60.
        self.dell_pad = 2*np.pi/(self.npad * reso)
        self.ell1d = 2*np.pi * np.fft.fftfreq(self.npad,d=reso)
        
        
    def expand_1d_to_2d(l,fl):
        lx = np.tile(self.ell1d, [1,self.npad])
        #should be npad x npad now
        ll = np.sqrt(lx*lx  + lx.T*lx.T)
        return np.interp(ll, l,fl,left=0,right=0)
    
    def camb_to_psd(self,spectra):
        psd = np.zeros([self.npad,self.npad])
        order = ['L','TT','EE','BB','TE']
        vec = range(5)
        oslice = order == spectra.upper()
        nfound = np.sum(oslice)
        if nfound is not 1:
            print('didnt find {}'.format(spectra))
            return np.zeros([self.npad,self.npad])
        else:
            ii = vec[order]
            Cl = self.cmb_spectra[ii,:]  #check that this is in Cl!
            ll = self.cmb_spectra[0,:]
            return expand_1d_to_2d(ll,Cl)
        
            
    def create_cmb_qu(self, Cov_4d = None):
        te_maps = self.create_cmb(Cov_4d=Cov_4d)

        tmap = np.squeeze(te_maps[:,:,0])
        qmap, umap = convert_e_to_qu(np.squeeze(te_maps[:,:,1]))
        return tmap, qmap, umap

    def create_cmb(self,Cov_4d = None):
        '''
        This is done in the TE basis.
        If you want IQU, a transform has to be applied afterwards.
        
        This will return a_kk's, The transfer function and beam should be applied afterwards (possibly after lensing)
        '''
        if Cov_4d is None:
            Cov_4d = np.zeros([self.npad,self.npad,2,2])
            theory_psd = self.camb_to_psd('TT')
            Cov_4d[:,:,0,0]=theory_psd
            theory_psd = self.camb_to_psd('EE')
            Cov_4d[:,:,1,1]=theory_psd
            theory_psd = self.camb_to_psd('TE')
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
        
        This will return a_kk's, The transfer function and beam should be applied afterwards (possibly after adding other terms).
        '''
        if Cov_4d is None:
            Cov_4d = np.zeros([self.npad,self.npad,self.nfreq,self.nfreq])
            k=0
            for i in range(self.nfreq):
                for j in range(i,self.nfreq):
                    theory_psd = self.fg_spectra_to_psd(i,j)
                    Cov_4d[:,:,i,j] = Cov_4d[:,:,j,i] = theory_psd
                    k+=1

        return create_correlated_power(Cov_4d)


    def create_noise(self,Cov_4d = None):
        '''
        Assumes have a set of PSDs
        An individual PSD is npad x npad in size
        For 3 observing freqs, have 6 PSDs --
        eg 90x90;, 90x150, 90x220, 150x150, 150x220, 220x220)
        ordering from self.pairs
        
        This will return a_kk's. This should be added to the signal terms after they have had filtering and beams applied. 
        (Presumably the PSD used here is from a filtered map and doesn't need additional transfer function due to filtering)
        '''
        if Cov_4d is None:
            Cov_4d = np.zeros([self.npad,self.npad,self.nfreq,self.nfreq])
            k=0
            for i in range(self.nfreq):
                for j in range(i,self.nfreq):
                    Cov_4d[:,:,i,j] = Cov_4d[:,:,j,i] = self.psds[k,:,:]
                    k+=1

        return create_correlated_power(Cov_4d)

    def lens_map(self, map, deflection_x_map, deflection_y_map, poly_deg = 5):
        '''
        Assume lensing is zero by the edge of the small map
        Could pass in larger map if desired.
        '''
        nx = deflection_x_map.shape[0]
        ny = deflection_x_map.shape[1]
        ynx = deflection_y_map.shape[0]
        yny = deflection_y_map.shape[1]
        mnx = maps.shape[0]
        mny = maps.shape[1]
        assert nx == ny == mnx == mny == yny == ynx  # only doing squares
        x = np.arange(0,nx) * (self.reso_arcmin * np.pi/180)
        x2d = np.tile(x,[1,nx])
        y2d = x2d.T
        
        mod_theta_x_grid = (x2d + deflection_x_map).flatten()
        mod_theta_y_grid = (x2d.T + deflection_y_map).flatten()

        lensed_map  = intrp.RectBivariateSpline( x, x, map, kx = poly_deg, ky = poly_deg).ev(mod_theta_y_grid, mod_theta_x_grid).reshape([nx,nx])
        return lensed_map

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



    def apply_beam(self,large_maps, index, map_in_fourier = True):
        '''Apply beams in Fourier space. Does not return to real space'''
        if not map_in_fourier:
            raise Exception("Haven't implemented real map in apply_beams")
        assert large.maps.shape[0] == self.npad
        assert large.maps.shape[1] == self.npad

        lbl = self.beams[0,:]
        bl = self.beams[index,:]
        bl2d = expand_1d_to_2d(lbl,bl)
        return large_maps * bl2d

    def convert_e_to_qu(fourier_e_largemap):
        '''
        Takes Fourier maps E (B assumed 0) and converts to Q/U (still in Fourier space)
        '''
        lx = np.tile(self.ell1d, [1,self.npad])

        angle = 2 * np.arctan2(lx, -1*lx.T)
        cosang = np.cos(angle)
        sinang = np.sin(angle)

        return cosang * fourier_e_largemap, sinang * fourier_e_largemap


    def create_sim_map(self):
        ft,fq,fu = self.create_cmb_qu()
        #get def angles
        #got to real space:
        t = np.fft.ifft2(ft).real
        q = np.fft.ifft2(fq).real
        u = np.fft.ifft2(fu).real
        #may only do the following on a portion of the large map...
        lens_t = self.lens_map(t, def_x, def_y)
        lens_q = self.lens_map(q, def_x, def_y)
        lens_u = self.lens_map(u, def_x, def_y)
        #and back to Fourier:
        f_lens_t = np.fft.fft2(lens_t)
        f_lens_q = np.fft.fft2(lens_q)
        f_lens_u = np.fft.fft2(lens_u)
        #this will be temp only
        fgs = self.create_fgs()
        # consider adding pol fgs, or continuing to assume zero
        #ufgs = self.create_fgs(pol=True)
        #qfgs = self.create_fgs(pol=True)

        ft_all = {}
        fq_all = {}
        fu_all = {}
        if i in range(self.nfreq):
            pdb.set_trace() # confirm if filtering depends on T/Q/U or not
            ft_all[self.freq_names[i]] = self.filter_map(self.apply_beam((f_lens_t + np.squeeze(fgs[:,:,0]) ) ,i),i)
            fq_all[self.freq_names[i]] = self.filter_map(self.apply_beam((f_lens_q ),i),i)
            fu_all[self.freq_names[i]] = self.filter_map(self.apply_beam((f_lens_u ) ,i),i)

        #could reorder this and above to reduce memory usage if useful
        noise = self.create_noise()
        out_t = {}
        out_q = {}
        out_u = {}
        if i in range(self.nfreq):
            pdb.set_trace() #almost certainly need to do something about how noise is defined
            tmp = ft_all[self.freq_names[i]] + tnoise

            pdb.set_trace() # need to confirm if desired output map is apodized or not
            out_t[freq_names[i]] = self.cut_map(np.fft.ifft2(tmp).real)  
            
            tmp = fq_all[self.freq_names[i]] + qnoise
            out_q[freq_names[i]] = self.cut_map(np.fft.ifft2(tmp).real) 

            tmp = fu_all[self.freq_names[i]] + unoise
            out_u[freq_names[i]] = self.cut_map(np.fft.ifft2(tmp).real) 
        
        return out_t, out_q, out_u