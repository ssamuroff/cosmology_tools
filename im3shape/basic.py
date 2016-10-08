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

    psfsize=(boxsize+opt._struct.contents.padding)*opt._struct.contents.upsampling
    upsampling=opt._struct.contents.upsampling

    large_boxsize = boxsize*6

    gal1 = get_model(size=size, type=model, g1=shear[0], g2=shear[1], flux=flux)
    pos1 = galsim.PositionD(0,0)

    # Add a tiny nominal PSF
    psf, psf_image = get_fixed_gaussian_psf(opt,psf_size, psf_ellipticity[0], psf_ellipticity[1], wcs)
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

def run(galaxy_image, psf_image, show=False, opt=None):
    if opt is None:
        opt = p3s.Options("/home/samuroff/shear_pipeline/end-to-end/end-to-end_code/config_files/im3shape/params_disc.ini")
    boxsize = galaxy_image.array.shape[0]
    wt=np.ones((boxsize,boxsize))
    trans = [transform(boxsize, opt)]
    result, best_img, images, weights = p3s.analyze_multiexposure([galaxy_image.array], [psf_image], [wt],trans, opt, ID=3000000, bands="?")

    if show:
        mosaic = np.vstack((galaxy_image.array,best_img))
        plt.imshow(mosaic,interpolation="none")
        plt.matshow(best_img)
        plt.colorbar()
        plt.matshow(galaxy_image.array)
        plt.colorbar()

    return result.get_params()

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
    def generate(self, dneigh=20, central_flux=1939.0, psf_size=3.6, neighbour_flux=2379.0, neighbour_size=3.0):
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
        
        self.psfe=np.linspace(-0.2,0.2,12)
        
        # Cycle over psf ellipticities
        for pe in self.psfe:
            g_theta["e1"]=[]
            g_theta["e2"]=[]
            ang=[]
            # For each psf ellipticity average over a ring of neighbour positions
            # For the moment we hardcode the other relevant parameters to their mean values from the simulation
            for i,t in enumerate(theta):
                x=dneigh*np.cos(t)
                y=dneigh*np.sin(t)
                gal,psf = setup_simple(boxsize=48,shear=(0,0.0), flux=central_flux, psf_ellipticity=(pe,0.0), psf_size=psf_size, neighbour_ellipticity=(0.0,0.0),neighbour=[x,y], neighbour_flux=neighbour_flux, neighbour_size=neighbour_size, wcs=wcs, opt=opt)
                res = run(gal,psf, opt=opt)
                g_theta["e1"].append(res.e1)
                gal,psf = setup_simple(boxsize=48,shear=(0,0.0), flux=central_flux, psf_ellipticity=(pe,0.0), psf_size=psf_size, neighbour_ellipticity=(0.0,0.0),neighbour=[np.inf,np.inf], neighbour_flux=neighbour_flux, neighbour_size=neighbour_size, wcs=wcs, opt=opt)
                res = run(gal,psf, opt=opt)
                g_theta["e2"].append(res.e1)
                ang.append(t)
            self.g["e1"].append(np.mean(g_theta["e1"]))
            self.g["e2"].append(np.mean(g_theta["e2"]))
    
        for k in self.g.keys(): self.g[k]=np.array(self.g[k])

    def show(self, xmin=-0.2,xmax=0.2,noneigh=True):

        fig = plt.figure()
        plt.plot(self.psfe,self.g["e1"]*1e4,color="purple",lw=2.5, ls="-",label="with neighbour")
        if noneigh:
            plt.plot(self.psfe,self.g["e2"]*1e4,color="steelblue",lw=2.5, ls=":",label="without neighbour")
        plt.axhline(0.0, color="k", linestyle="--", lw=2.5)
        plt.axvline(0.0, color="k", linestyle="--", lw=2.5)
        plt.xlim(xmin,xmax)

        plt.xlabel(r"PSF ellipticity $ e_1 ^{psf}$")
        plt.ylabel(r"mean ellipticity $<e_1>_{\theta} \times 10^4$")
        if noneigh:
            plt.legend(loc="upper left")
        plt.show()

    def alpha(self, component=1):
        sel = abs(self.psfe)<0.1
        linear_e = self.g["e%d"%component][sel]
        linear_epsf = self.psfe[sel]
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


gs = gridspec.GridSpec(3,1)
ax1 = plt.subplot(gs[:2,:])
ax2 = plt.subplot(gs[2:,:])
ax1.plot(x,alpha, "-",lw=2.5,color="purple")
ax2.hist(sim.res["mean_psf_fwhm"], bins=np.linspace(x[0],x[-1],50), histtype="step", color="steelblue", lw=2.5, normed=1)
plt.subplots_adjust(wspace=0, hspace=0)
t=[1,2,3,4,5]
ax1.xticks([],[])
ax2.xticks(t,["%1.1f"%i for i in t])
ax2.yticks([],[])
ax2.add_ylabel(r"PSF leakage $\alpha _1$")
ax2.add_xlabel(r"PSF FWHM / pixels")
plt.show()







