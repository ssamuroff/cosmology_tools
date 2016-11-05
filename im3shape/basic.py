import numpy as np
import galsim
from py3shape import structs
import math
import pylab as plt
import fitsio
import glob, argparse, os
import tools.shapes as s
import py3shape as p3s
from py3shape import utils

sersic_indices={"disc":1, "bulge":4}

labels={"psf_size":("mean_psf_fwhm", r"PSF FWHM / pixels"), "neighbour_flux":("mean_flux",r"Neighbour Flux $f_{neigh}$"), "central_flux":("mean_flux",r"Central Flux $f_{cent}$"), "neighbour_size":("radius",r"Neighbour Size $R_{neigh}$/ pixels")}

def get_model(type="disc", size=1.0, g1=0, g2=0, flux=None):
    if type is "disc":
        gal = galsim.Sersic(n=1, half_light_radius=size)
    if type is "bulge":
        gal = galsim.Sersic(n=4, half_light_radius=size)
    if type is "gaussian":
        gal = galsim.Gaussian(fwhm=size)

    shear = galsim.Shear(g1=g1,g2=g2)
    gal = gal.shear(shear)

    if flux is not None:
        gal = gal.withFlux(flux)

    return gal


def get_fixed_gaussian_psf(options, psf_size, psf_e1, psf_e2, wcs=None):
    "This is slow - for test runs only."
    import galsim
    psf_box=(options.stamp_size+options.padding)*options.upsampling
    print psf_box

    #Get the localized WCS information
    if wcs is None:
        wcs_path = "/share/des/disc2/y1/OPS/coadd/20141118000051_DES0014-4414/coadd/DES0014-4414_r.fits.fz"
        wcs = galsim.FitsWCS(wcs_path)

    orig_col = 1000
    orig_row = 1000
    image_pos = galsim.PositionD(orig_col,orig_row)
    local_wcs = wcs.local(image_pos)
    local_wcs._dudx = 1.0/options.upsampling
    local_wcs._dudy = 0.0/options.upsampling
    local_wcs._dvdx = 0.0/options.upsampling
    local_wcs._dvdy = 1.0/options.upsampling
    A=np.array([[local_wcs._dudx,local_wcs._dudy],[local_wcs._dvdx, local_wcs._dvdy]])
    local_wcs._det =np.linalg.det(A)

    #Generate the PSF in sky coordinates
    pix = wcs.toWorld(galsim.Pixel(0.27), image_pos=image_pos)
    psf_sky = galsim.Convolve([galsim.Gaussian(fwhm=psf_size), pix])
    psf_sky = galsim.Gaussian(fwhm=psf_size)
    pshear = galsim.Shear(g1=psf_e1,g2=psf_e2)
    psf_sky = psf_sky.shear(pshear)

    psf_stamp = psf_sky.drawImage(wcs=local_wcs, nx=psf_box, ny=psf_box, method='no_pixel')
    psf_stamp_array = psf_stamp.array.copy()
    assert psf_stamp_array.shape == (psf_box, psf_box)
    return psf_sky, psf_stamp.array.copy()

def setup_simple(boxsize=32, model="disc", size=2, flux=None, shear=(0,0), neighbour=[np.inf,np.inf], neighbour_size=2.0, neighbour_flux=0.13, neighbour_ellipticity=(0,0),psf_ellipticity=(0,0), psf_size=0.1, wcs=None, opt=None):
    """Basic function to construct a test galaxy for exploring shape biases."""

    if opt is None:
        opt=p3s.Options("/home/samuroff/shear_pipeline/end-to-end/end-to-end_code/config_files/im3shape/params_disc.ini")

    opt.stamp_size= boxsize

    psfsize=(boxsize+opt._struct.contents.padding)*opt._struct.contents.upsampling
    upsampling=opt._struct.contents.upsampling

    print psfsize

    large_boxsize = boxsize*6

    gal1 = get_model(size=size, type=model, g1=shear[0], g2=shear[1], flux=flux)
    pos1 = galsim.PositionD(0,0)

    # Add a tiny nominal PSF
    psf, psf_image = get_fixed_gaussian_psf(opt,psf_size, psf_ellipticity[0], psf_ellipticity[1], wcs)
    print psf_image.shape
    #import pdb ; pdb.set_trace()
    #galsim.Gaussian(fwhm=psf_size*0.27)

    if wcs is None:
        wcs_path = "/share/des/disc2/y1/OPS/coadd/20141118000051_DES0014-4414/coadd/DES0014-4414_r.fits.fz"
        wcs = galsim.FitsWCS(wcs_path)
    orig_col = 1000
    orig_row = 1000
    image_pos = galsim.PositionD(orig_col,orig_row)
    local_wcs = wcs.local(image_pos)

    local_wcs._dudx = 1.0
    local_wcs._dudy = 0.0
    local_wcs._dvdx = 0.0
    local_wcs._dvdy = 1.0
    A=np.array([[local_wcs._dudx,local_wcs._dudy],[local_wcs._dvdx, local_wcs._dvdy]])
    local_wcs._det =np.linalg.det(A)

    pix = wcs.toWorld(galsim.Pixel(0.27), image_pos=image_pos)
    
    obj1= galsim.Convolve([gal1,psf,pix])

    im1 = obj1.drawImage(wcs=local_wcs, nx=large_boxsize, ny=large_boxsize, method='no_pixel')
    im3 = obj1.drawImage(wcs=local_wcs, nx=large_boxsize, ny=large_boxsize, method='no_pixel')

    # And a neighbouring object if specified
    if np.isfinite( np.array(neighbour)).all():
        print "drawing neighbour at (%d, %d)"%(int(neighbour[0]), int(neighbour[1]))
        gal2=galsim.Sersic(n=1,half_light_radius=neighbour_size)
        gal2=gal2.withFlux(neighbour_flux)
        nshear=galsim.Shear(g1=neighbour_ellipticity[0],g2=neighbour_ellipticity[1])
        gal2=gal2.shear(nshear)

        obj2= galsim.Convolve([gal2,psf,pix])

        pos2 = galsim.PositionD(x=int(neighbour[0]), y=int(neighbour[1]))

        im2 = obj2.drawImage(wcs=local_wcs, nx=large_boxsize, ny=large_boxsize, method='no_pixel', offset=pos2)

        im3+=im2.array

    trim = (large_boxsize-boxsize)/2
    lim = galsim.BoundsI(xmin=trim, xmax=large_boxsize-trim-1, ymin=trim, ymax=large_boxsize-trim-1)

    return im3[lim], psf_image

def run(galaxy_image, psf_image, save=None, show=False, weights=[[None]], opt=None, ncols=1, col=1, return_cat=False, flat_background=0, noise=0, check=False):
    if opt is None:
        opt = p3s.Options("/home/samuroff/shear_pipeline/end-to-end/end-to-end_code/config_files/im3shape/params_disc.ini")
    boxsize = galaxy_image.array.shape[0]
    opt.stamp_size = boxsize
    if weights[0][0] is None:
        wt=np.ones((boxsize,boxsize))
    else:
        wt = weights
    trans = [transform(boxsize, opt)]
    if noise==0:
        noise= np.zeros_like(wt)
    else:
        noise=noise*(np.random.rand(wt.size).reshape(wt.shape[0],wt.shape[1])-0.5)
    result, best_img, images, weights = p3s.analyze_multiexposure([galaxy_image.array+flat_background+noise], [psf_image], [wt],trans, opt, ID=3000000, bands="?")

    if show:
        cmap = plt.get_cmap("jet")
        mosaic = np.vstack((galaxy_image.array,best_img))

        plt.subplot(2,ncols, col)
        plt.imshow(galaxy_image.array, interpolation="none", cmap=cmap)
        plt.subplot(2,ncols,ncols+col)
        plt.imshow(best_img, interpolation="none", cmap=cmap)
        plt.subplots_adjust(hspace=0)

    if save is not None:
        np.savetxt(save+"/model_tm.txt", best_img)
        np.savetxt(save+"/image_tm.txt", galaxy_image.array)
        np.savetxt(save+"/weights_tm.txt", wt)

    if check:
        sel = check_cuts(result)
        print "Object passes info_flag cuts:", sel
        return result.get_params(), result, check


    if return_cat:
        return result.get_params(), result
    else:
        return result.get_params()

def check_cuts(results):
    param = results.get_params()

    ra_as = param.ra_as
    dec_as = param.dec_as
    R = np.sqrt(ra_as*ra_as + dec_as*dec_as)
    snr = results.snr

    e = [param.e1, param.e2]
    radius = param.radius
    min_residuals = results.min_residuals
    max_residuals = results.max_residuals

    print "passes SNR = %f > 10 : "%snr, snr > 10.
    print "passes SNR = %f < 10000 : "%snr, snr < 10000.
    print "passes R = %f < 0.6 arcsec : "%R, R < 0.6/0.27
    print "passes radius = %f < 5.0 : "%radius, radius < 5.0
    print "passes radius = %f> 0.1 : "%radius, radius > 0.1
    print "passes -1 < e1 = %f < 1 : "%e[0], -1.0 < e[0] < 1.0
    print "passes -1 < e2 = %f < 1 : "%e[1], -1.0 < e[1] < 1.0
    print "passes min_residuals = %f > -0.2 : "%min_residuals, min_residuals > -0.2
    print "passes max_residuals = %f < 0.2 : "%max_residuals, max_residuals < 0.2

    return (snr > 10.) & (snr < 10000.) & (R < 0.6/0.27) & (radius < 5.0) & (radius > 0.1) & (min_residuals > -0.2) & (max_residuals < 0.2) & (-1.0 < e[0] < 1.0) & (-1.0 < e[1] < 1.0)

    #info_cuts =["(%s['fails_unmasked_flux_frac']==0)", "(%s['snr']>10)", "(%s['snr']<10000)", "(%s['mean_rgpp_rp']>1.1)", "(%s['mean_rgpp_rp']<3.5)", "(%s['radius']<5)", "(%s['radius']>0.1)", "((%s['ra_as']**2+%s['dec_as']**2)**0.5<1.0)", "(%s['chi2_pixel']>0.5)", "(%s['chi2_pixel']<1.5)", "(%s['min_residuals']>-0.2)", "(%s['max_residuals']<0.2)", "(%s['mean_psf_fwhm']<7.0)", "(%s['mean_psf_fwhm']>0.0)", "(%s['error_flag']==0)"]



def transform(boxsize, opt=None):
    """Dummy im3shape transform, setting the MCMC chain started in the centre of the stamp."""

    if opt is None:
        options = p3s.Options("/home/samuroff/shear_pipeline/end-to-end/end-to-end_code/config_files/im3shape/params_disc.ini")
        
    transform = structs.Transform_struct()
    transform.ra0 = 0.0  # these are not used in this case
    transform.dec0 = 0.0
    transform.cosdec0 = 0.0 
    transform.x0 = boxsize/2
    transform.y0 = boxsize/2
    print "starting chain with centroid position (%d %d)"% (transform.x0 ,transform.y0 )
    dudcol=1.0
    dudrow=0.0
    dvdcol=0.0
    dvdrow=1.0
    B = np.array([ [dudcol,dudrow],[dvdcol,dvdrow]])
    A = np.linalg.inv(B)
    transform.A[0][0] = A[0,0]
    transform.A[0][1] = A[0,1]
    transform.A[1][0] = A[1,0]
    transform.A[1][1] = A[1,1]
    #print 'transform det = ', np.linalg.det(B)
    #print 'A = ', A
    return transform

class toy_model:
    def generate(self, dneigh=20, central_flux=1912.0, psf_size=3.7, neighbour_flux=1912.0, neighbour_size=3.2, vary="psf_e"):
        # PSF leakage toy model
        #-----------------------------------------------------------------------
        g_theta={}
        theta= np.linspace(0,2,50)*np.pi
        self.g={}
        self.g["e1"]=[]
        self.g["e2"]=[]
        ang=[]
    
        # Setup a dummy wcs
        # We load it here to avoid too much unncessary io
        # In fact it should provide a skeleton only, since we overwrite the Jacobian with an identity matrix
        wcs_path = "/share/des/disc2/y1/OPS/coadd/20141118000051_DES0014-4414/coadd/DES0014-4414_r.fits.fz"
        orig_col = 1000
        orig_row = 1000
        image_pos = galsim.PositionD(orig_col,orig_row)
        wcs = galsim.FitsWCS(wcs_path)
        opt = p3s.Options("/home/samuroff/shear_pipeline/end-to-end/end-to-end_code/config_files/im3shape/params_disc.ini")
        
        binning = {"psf_e1": np.linspace(-0.2,0.2,12), "cent_e1": [-0.65,0.65]}
        self.binning = binning[vary]
        
        # Cycle over psf ellipticities
        for e in self.binning:
            if vary=="psf_e1":
                pe = e
                ce = 0
            elif vary=="cent_e1":
                pe = 0
                ce = e
            g_theta["e1"]=[]
            g_theta["e2"]=[]
            ang=[]
            # For each psf ellipticity average over a ring of neighbour positions
            # For the moment we hardcode the other relevant parameters to their mean values from the simulation
            for i,t in enumerate(theta):
                x=dneigh*np.cos(t)
                y=dneigh*np.sin(t)
                gal,psf = setup_simple(boxsize=48,shear=(ce,0.0), flux=central_flux, psf_ellipticity=(pe,0.0), psf_size=psf_size, neighbour_ellipticity=(0.0,0.0),neighbour=[x,y], neighbour_flux=neighbour_flux, neighbour_size=neighbour_size, wcs=wcs, opt=opt)
                res = run(gal,psf, opt=opt)
                g_theta["e1"].append(res.e1)
                gal,psf = setup_simple(boxsize=48,shear=(ce,0.0), flux=central_flux, psf_ellipticity=(pe,0.0), psf_size=psf_size, neighbour_ellipticity=(0.0,0.0),neighbour=[np.inf,np.inf], neighbour_flux=neighbour_flux, neighbour_size=neighbour_size, wcs=wcs, opt=opt)
                res = run(gal,psf, opt=opt)
                g_theta["e2"].append(res.e1)
                ang.append(t)
            self.g["e1"].append(np.mean(g_theta["e1"]))
            self.g["e2"].append(np.mean(g_theta["e2"]))
    
        for k in self.g.keys(): self.g[k]=np.array(self.g[k])

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


class toy_model2:
    def __init__(self, komplexitet=1):
        print "Toy model to illustrate neighbour bias."
        print "Complexity level: %d"%komplexitet
        self.complexity=komplexitet

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


    def analyse1(self, central_ellipticity, dneigh=20, central_flux=1912.0, psf_size=3.7, neighbour_flux=1912.0, neighbour_size=3.2):
        # Level 1 - loop over SNR
        #-----------------------------------------------------------------------
        
        self.m={}
        self.m[1]=[]
        self.m[2]=[]
        self.snr=[]

        self.dneigh=20
        self.central_flux=1912.0
        self.psf_size=3.7
        self.neighbour_flux=1912.0
        self.neighbour_size=3.2

        print "Toy model has a central galaxy with shape ", central_ellipticity
        self.e1 = central_ellipticity[0]
        self.e2 = central_ellipticity[1]
    
        # Setup a dummy wcs
        # We load it here to avoid too much unncessary io
        # In fact it should provide a skeleton only, since we overwrite the Jacobian with an identity matrix
        wcs_path = "/share/des/disc2/y1/OPS/coadd/20141118000051_DES0014-4414/coadd/DES0014-4414_r.fits.fz"
        orig_col = 1000
        orig_row = 1000
        image_pos = galsim.PositionD(orig_col,orig_row)
        self.wcs = galsim.FitsWCS(wcs_path)
        self.opt = p3s.Options("/home/samuroff/shear_pipeline/end-to-end/end-to-end_code/config_files/im3shape/params_disc.ini")
        
        self.binning = np.logspace(1,3.5,12)
        
        # Cycle over the central flux, as a proxy for SNR
        for i, proxy in enumerate(self.binning):
            print "Level 1 iteration: %d %f"%(i,proxy)
            self.do_position_loop(proxy)
        
        print "Done all loops"

    def do_position_loop(self, proxy):
        """Loop over all possible angular positions for this neighbour realisation."""
        ang=[]
        m_theta={}
        theta= np.linspace(0,2,10)*np.pi
        snr = []

        for i,t in enumerate(theta):
            print "  --  Level 2 iteration: %d (theta=%f)"%(i,t)

            m_theta[1]=[]
            m_theta[2]=[]

            x = self.dneigh*np.cos(t)
            y = self.dneigh*np.sin(t)

            sn, m1,m2 = self.simple_m(proxy, x, y)
            print "At SNR = %4.1f got m = (%2.2f,%2.2f)"%(sn,m1,m2)
            m_theta[1].append(m1)
            m_theta[2].append(m2)
            ang.append(t)
            snr.append(sn)

        # Compute the mean over all positions
        self.m[1].append(np.mean(m_theta[1]))
        self.m[2].append(np.mean(m_theta[2]))
        self.snr.append(np.mean(snr))

    def simple_m(self, proxy, x, y):

        y_g = []
        snr = []
        shear = [-0.02,0.02]
        for g in shear:
            gal, psf = setup_simple(boxsize=48, shear=(self.e1+g,self.e2), flux=proxy, psf_ellipticity=(0.0,0.0), psf_size=self.psf_size, neighbour_ellipticity=(0.0,0.0), neighbour=[x,y], neighbour_flux=self.neighbour_flux, neighbour_size=self.neighbour_size, wcs=self.wcs, opt=self.opt)
            res, cat = run(gal, psf, opt=self.opt, return_cat=True, show=False)
            snr.append(cat.snr)
            y_g.append(res.e1 - self.e1 - g)

        m1 = (y_g[1]-y_g[0])/(shear[1]-shear[0])
        snr = np.mean(snr)

        y_g = []
        for g in shear: 
            gal, psf = setup_simple(boxsize=48, shear=(self.e1,self.e2+g), flux=proxy, psf_ellipticity=(0.0,0.0), psf_size=self.psf_size, neighbour_ellipticity=(0.0,0.0), neighbour=[x,y], neighbour_flux=self.neighbour_flux, neighbour_size=self.neighbour_size, wcs=self.wcs, opt=self.opt)
            res, cat = run(gal, psf, opt=self.opt, return_cat=True)
            y_g.append(res.e2 - self.e2 -  g)

        m2 = (y_g[1]-y_g[0])/(shear[1]-shear[0])

        return snr, m1, m2

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












def plot_alpha(filename, sim, show_mean=True, points=False, ncat=None):

    from matplotlib import gridspec

    param_name = os.path.basename(filename).replace("_toymodel_data2.txt","")

    x,alpha=np.loadtxt(filename).T
    dx=(x[1]-x[0])

    if param_name!="dneigh":
        label = labels[param_name][1]
        sim_distn = sim.res[ labels[param_name][0] ]
    else:
        label = r"Neighbour distance $d_{neigh} / R_{stamp}$"
        dr = (sim.truth["ra"]-ncat.truth["ra"])**2
        dr+= (sim.truth["dec"]-ncat.truth["dec"])**2
        dr = np.sqrt(dr)
        sim_distn = dr*60*60/0.27 / sim.res["stamp_size"]
        x/=48.0

    gs = gridspec.GridSpec(3,1)

    ax1 = plt.subplot(gs[:2,:])
    ax2 = plt.subplot(gs[2:,:])
    ax1.plot(x,alpha, "-",lw=2.5,color="purple")
    if points:
        ax1.plot(x,alpha, "o",lw=2.5,color="purple")
    ax2.hist( sim_distn, bins=np.linspace(x[0],x[-1],45), histtype="step", color="steelblue", lw=2.5, normed=1)

    plt.subplots_adjust(wspace=0, hspace=0)
 

    t=np.linspace(x[0],x[-1],6)
    ax1.set_xticks([],[])
    ax2.set_xticks(t,["%1.1f"%i for i in t])
    ax2.set_yticks([],[])
    ax1.set_ylabel(r"PSF leakage $\alpha _1$")
    ax2.set_xlabel(label)

    ax1.set_xlim(x[0],x[-1])
    ax2.set_xlim(x[0],x[-1])

    if show_mean:
        ax1.axvline(sim_distn.mean(), color="k", ls="--")
        ax2.axvline(sim_distn.mean(), color="k", ls="--")

    plt.show()






