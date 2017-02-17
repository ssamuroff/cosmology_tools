
import numpy as np
import galsim
import fitsio
import glob, argparse, os
import tools.diagnostics as di
import tools.plots as pl
import fitsio as fi
import pylab as plt
import os, yaml, argparse
import tools.shapes as s


description = 'Noise free toy model'
parser = argparse.ArgumentParser(description=description, add_help=False)
parser.add_argument('--nrealisations', type=int, action='store', default=2000)
parser.add_argument('--output', type=str, default="/home/samuroff/hoopoe_paper/toy_model_data")
parser.add_argument('--fatal_errors', action='store_true')
parser.add_argument('--mpi', action='store_true')
args = parser.parse_args()

# Mean neighbour distance (pixels) : 39.79


print "Will calculate m samples"
print "-------------------------------------------"


try:
	print "Setting up MPI"
	import mpi4py.MPI
	rank = mpi4py.MPI.COMM_WORLD.Get_rank()
	size = mpi4py.MPI.COMM_WORLD.Get_size()

except:
	rank=0
	size=1

config = yaml.load(open("/home/samuroff/calibration_config.yaml"))
hoopoe = s.shapecat(res=config["hoopoe_dir"], truth=config["hoopoe_dir"] )
hoopoe.res = fi.FITS(config["hoopoe_dir"])["i3s"][:args.nrealisations*4]
hoopoe.truth = fi.FITS(config["hoopoe_dir"])["truth"][:args.nrealisations*4]
sel=((hoopoe.res["snr"] > 12) & (hoopoe.res["snr"] < 200) & (hoopoe.res["mean_rgpp_rp"] > 1.13) & (hoopoe.res["mean_rgpp_rp"] < 3.0))
hoopoe.res=hoopoe.res[sel]
hoopoe.truth=hoopoe.truth[sel]

ncat=fi.FITS("/home/samuroff/neighbours.fits")[1].read()[:args.nrealisations*4]

f=np.sqrt(-2*np.log(0.5))                    # HLR - sigma conversion
fac=2.35                                     # sigma - FWHM conversion


#ets=np.sqrt( (hoopoe.truth["intrinsic_e1"]+hoopoe.truth["true_g1"])**2 + (hoopoe.truth["intrinsic_e2"]+hoopoe.truth["true_g2"])**2 )
et=np.sqrt( (hoopoe.truth["intrinsic_e1"])**2 + (hoopoe.truth["intrinsic_e2"])**2 )
R=hoopoe.truth["hlr"]/f                      #Gaussian sigma galaxy, arcsec
rp=hoopoe.res["mean_hsm_psf_sigma"]*0.27/3   # sigma PSF, factor of 3 for the upsampling, arcsec
#rp=rp.mean()

inp=np.empty(hoopoe.res.size, dtype=[("hlr", float),("psf_size", float),("flux", float), ("neighbour_flux", float), ("neighbour_hlr", float), ("nearest_neighbour_pixel_dist", float)])
inp["nearest_neighbour_pixel_dist"]=hoopoe.truth["nearest_neighbour_pixel_dist"]
f=np.sqrt(-2*np.log(0.5))                    # HLR - sigma conversion
fac=2.35                                     # sigma - FWHM conversion


R=hoopoe.truth["hlr"]/f                      #Gaussian sigma galaxy, arcsec
rp=hoopoe.res["mean_hsm_psf_sigma"]*0.27/3   # sigma PSF, factor of 3 for the upsampling, arcsec
inp["psf_size"]=rp*fac/0.27
inp["hlr"]=R*fac/0.27
inp["flux"]=hoopoe.truth["flux"]
inp["neighbour_flux"]=ncat["flux"]
inp["neighbour_hlr"]=ncat["hlr"]

from tools.im3shape import mcmc_toy_model as mc
model=mc.mcmc_toy_model() 
model.run(inp, filename="%s/mc_toy_model-results-%d.txt"%(args.output,rank), size=size, rank=rank)