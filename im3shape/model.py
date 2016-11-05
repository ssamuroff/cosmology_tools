import galsim, os
import numpy as np
import pylab as plt
from tools.im3shape import basic as i3s
import tools.diagnostics as di
import py3shape as p3s

class toy_model2:
    def __init__(self, komplexitet=1, random_seed=90000):
        print "Toy model to illustrate neighbour bias."
        print "Complexity level: %d"%komplexitet
        self.complexity=komplexitet
        self.random_seed = random_seed

    def explain(self):
        print explanations[self.complexity]

    def set_complexity(self, komplexitet):
        print "Changed complexity level from %d to %d"%(self.complexity, komplexitet)
        self.complexity=komplexitet

    def show_params(self):
        print "Toy model has input parameters:"
        print "Central flux : %3.2f"%self.central_flux
        print "PSF size : %3.2f"%self.psf_size
        print "Neighbour distance : %3.2f"%self.dneigh
        print "Neighbour flux : %3.2f"%self.neighbour_flux
        print "Neighbour size : %3.2f"%self.neighbour_size


    def analyse1(self, central_data, neighbour_data, nreal=1000): #, central_ellipticity, dneigh=20, central_flux=1912.0, psf_size=3.7, neighbour_flux=1912.0, neighbour_size=3.2, nrealisations=1000):
        # Level 1 - loop over SNR
        #-----------------------------------------------------------------------
        
        # Initialise the random number generator
        np.random.seed(self.random_seed)

        self.object_mask = np.loadtxt("mask_template.txt")

        m={}
        m[1]=[]
        m[2]=[]
        c={}
        c[1]=[]
        c[2]=[]
        snr=[]

    
        # Setup a dummy wcs
        # We load it here to avoid too much unncessary io
        # In fact it should provide a skeleton only, since we overwrite the Jacobian with an identity matrix
        wcs_path = "/share/des/disc2/y1/OPS/coadd/20141118000051_DES0014-4414/coadd/DES0014-4414_r.fits.fz"
        orig_col = 1000
        orig_row = 1000
        image_pos = galsim.PositionD(orig_col,orig_row)
        self.wcs = galsim.FitsWCS(wcs_path)
        self.opt = p3s.Options("/home/samuroff/shear_pipeline/end-to-end/end-to-end_code/config_files/im3shape/params_disc.ini")
        
        self.binning = np.logspace(1,2.8,12)
        upper = self.binning[1:]
        lower = self.binning[:-1]
        
        # Cycle over the central flux, as a proxy for SNR
        for i, limits in enumerate(zip(lower, upper)):
            snr_min = limits[0]
            snr_max = limits[1]
            print "Level 1 iteration: %d SNR = %3.3f - %3.3f"%(i+1,snr_min, snr_max)
            # Selection function for this bin
            sel = (central_data.res["snr"]>snr_min) & (central_data.res["snr"]<snr_max)

            print "Will do %d realisations for this bin:"%nreal

            De1=[]
            De2=[]
            sn=[]
            g1=[]
            g2=[]
            for j in xrange(nreal):
                print "%d/%d"%(j+1,nreal)
                self.generate_central_realisation(central_data, sel)
                self.generate_neighbour_realisation(neighbour_data, central_data, sel)

                # Draw this neighbour realisation repeatedly on a ring of angular positions
                snr, e1, e2 = self.do_position_loop()

                sn.append(snr)
                De1.append(e1)
                De2.append(e2)
                g1.append(self.g[0])
                g2.append(self.g[1])

            data = np.zeros(len(g1), dtype=[("e1", float), ("e2", float), ("true_g1", float), ("true_g2", float)])
            data["e1"], data["e2"] = np.array(De1), np.array(De2) 
            data["true_g1"], data["true_g2"] = np.array(g1), np.array(g2) 

            bias = di.get_bias(data, nbins=5, names=["m11","m22","c11","c22"], binning="equal_number", silent=True)
            print bias

            m[1].append(bias["m11"][0])
            m[2].append(bias["m22"][0])
            c[1].append(bias["c11"][0])
            c[2].append(bias["c22"][0])
            
        import pdb ; pdb.set_trace()
        
        print "Done all loops"

    def generate_neighbour_realisation(self, neighbour_data, central_data, selection):
        print "Creating neighbour realisation ",
        i = np.random.randint(neighbour_data.truth[selection].size)
        print "random index: %d"%i

        self.neighbour_distance = neighbour_data.truth["Rn"][selection][i]
        self.neighbour_radius = neighbour_data.truth["hlr"][selection][i]
        self.neighbour_flux = neighbour_data.truth["flux"][selection][i]
        self.neighbour_e = (neighbour_data.truth["intrinsic_e1"][selection][i], neighbour_data.truth["intrinsic_e2"][selection][i])

        # Randomise the shear
        self.g = (np.random.rand(2)*0.1 - 0.05)

    def generate_central_realisation(self, central_data, selection):
        print "Creating central galaxy realisation."
        if self.complexity==1:
            print "Will use the mean values in this SNR bin."
            self.central_flux = central_data.truth["flux"][selection].mean()
            self.psf_e = (0,0)
            self.psf_size = central_data.res["mean_psf_fwhm"][selection].mean()
            self.central_e = (0,0)
            self.central_radius = central_data.res["radius"][selection].mean()

        elif self.complexity>1:
            print "Will use draw a random set of values in the relevant properties."
            i = np.random.randint(central_data.truth[selection].size)
            print "random index: %d"%i
            self.central_flux = central_data.truth["flux"][selection][i]
            self.psf_e = (0,0)
            self.psf_size = central_data.res["mean_psf_fwhm"][selection][i]
            self.central_e = (0,0)
            self.central_radius = central_data.res["radius"][selection][i]


    def do_position_loop(self):
        """Loop over all possible angular positions for this neighbour realisation."""
        # Some empty datavectors
        De1 = []
        De2 = []
        snr = []

        theta= np.linspace(0,2,10)*np.pi

        for i,t in enumerate(theta):
            print "  --  Level 2 iteration: %d (theta=%f)"%(i,t)

            x = self.neighbour_distance*np.cos(t)
            y = self.neighbour_distance*np.sin(t)

            gal, psf = i3s.setup_simple(boxsize=32, shear=self.g, flux=self.central_flux, psf_ellipticity=self.psf_e, psf_size=self.psf_size, neighbour_ellipticity=self.neighbour_e, neighbour=[x,y], neighbour_flux=self.neighbour_flux, neighbour_size=self.neighbour_radius, wcs=self.wcs, opt=self.opt)
            if self.neighbour_distance< 23:
                wts = self.generate_mock_weightmap(gal.array, x, y)
            else:
                wts = np.ones_like(gal.array)
            res, cat = i3s.run(gal, psf, weights=wts, opt=self.opt, return_cat=True, show=False)

            De1.append(res.e1)
            De2.append(res.e2)
            snr.append(cat.snr)

        # Compute the mean over all positions
        return np.mean(snr), np.mean(De1), np.mean(De2)

    def generate_mock_weightmap(self, image, x, y):
        """Use construct a mock seg map using a mask template.
           Then run uberseg to get a weight map for the stamp. """

        #Create a padded array for the seg map 
        newseg = np.zeros((image.shape[0]+20,image.shape[1]+20))
        newseg = np.zeros_like(newseg)

        n0 = newseg.shape[0]/2
        # Fill in the cental object and neighbour masks
        newseg[n0-self.object_mask.shape[0]/2:n0+self.object_mask.shape[0]/2,n0-self.object_mask.shape[1]/2:n0+self.object_mask.shape[1]/2]+=self.object_mask
        newseg[n0+int(y)-self.object_mask.shape[0]/2:n0+int(y)+self.object_mask.shape[0]/2,n0+int(x)-self.object_mask.shape[1]/2:n0+int(x)+self.object_mask.shape[1]/2]+=self.object_mask*2

        # Then trim the stamp
        newseg=newseg[10:-10,10:-10]

        return get_cweight_cutout_nearest(image,newseg)

    

    def show(self, xmin=-0.2,xmax=0.2,noneigh=True):

        fig = plt.figure()
        plt.plot(self.binning,self.g["e1"]*1e4,color="purple",lw=2.5, ls="-",label="with neighbour")
        if noneigh:
            plt.plot(self.binning,self.g["e2"]*1e4,color="steelblue",lw=2.5, ls=":",label="without neighbour")
        plt.axhline(0.0, color="k", linestyle="--", lw=2.5)
        plt.axvline(0.0, color="k", linestyle="--", lw=2.5)
        plt.xlim(xmin,xmax)

        plt.xlabel(r"PSF ellipticity $ e_1 ^{psf}$")
        plt.ylabel(r"mean ellipticity $<e_1>_{\theta} \times 10^4$")
        if noneigh:
            plt.legend(loc="upper left")
        plt.show()

    def alpha(self, component=1):
        sel = abs(self.binning)<0.1
        linear_e = self.g["e%d"%component][sel]
        linear_epsf = self.binning[sel]
        return (linear_e[-1]-linear_e[0])/(linear_epsf[-1]-linear_epsf[0])

    def sample_alpha(self, param_name, samples):
        alpha = []
        print "Will sample alpha at points in %s"%param_name
        for q in samples:
            print q,
            exec "self.generate(%s=%f)"%(param_name,q)
            alpha.append(self.alpha())
            print alpha[-1]

        return samples, alpha


def get_cweight_cutout_nearest(im,seg):

    weight = np.ones_like(im)
    obj_inds = np.where(seg != 0)
    # the seg map holds the sextractor number, 1 offset
    object_number = seg[seg.shape[0]/2,seg.shape[1]/2]
    
    # Then loop through pixels in seg map, check which obj ind it is closest
    # to.  If the closest obj ind does not correspond to the target, set this
            # pixel in the weight map to zero.
    for i,row in enumerate(seg):
        for j, element in enumerate(row):
            obj_dists = (i-obj_inds[0])**2 + (j-obj_inds[1])**2
            ind_min=np.argmin(obj_dists)

            segval = seg[obj_inds[0][ind_min],obj_inds[1][ind_min]]
            if segval != object_number:
                weight[i,j] = 0.

    return weight


exp1 = """Complexity level 1 --------------------------------------
  -- Drawing from the observed properies in each bin to 
     generate neighbour realisations.
  -- The properties of the central galaxy are fixed to 
     the means of the relvant distributions in that bin
  -- Using im3shape in disc mode
Assumptions:
  -- Round Gaussian PSFs 
  -- Round Central Galaxy
  -- No pixel noise
  -- Uniform random shear g1,2 = [-0.05,0.05]

  -----------------------------------------------------------------"""

explanations={1:exp1}