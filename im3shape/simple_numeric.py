import numpy as np
import galsim
from py3shape import structs
import math
import pylab as plt
plt.switch_backend("agg")
import fitsio
import glob, argparse, os
import tools.shapes as s
import py3shape as p3s
from py3shape import utils
from tools.im3shape import basic as i3s

Rc=2.1
Rp=1.3
Rn=2.9
fc=1650
fn=473

Rc_med=2.1
fc_med=945
fn_med=475
Rn_med=1.47

m=s.meds_wrapper("/share/des/disc8/cambridge/meds/DES2111+0043-r-sim-ohioA6-meds-y1a1-beta.fits.fz")

def make_sersic(e1,e2):
    gal = i3s.setup_simple(boxsize=32,shear=(e1,e2), psf_size=Rp,  size=Rc_med, neighbour_ellipticity=(0.0,0.0), neighbour_flux=fn_med, flux=fc_med, neighbour_size=Rn_med, neighbour=[np.inf,0], opt=m.options)
    return gal[0].array

g,p=i3s.setup_simple(boxsize=32,shear=(0.05,0.2), psf_size=Rp,  size=Rc_med, neighbour_ellipticity=(0.0,0.0), neighbour_flux=fn_med, flux=fc_med, neighbour_size=Rn_med, neighbour=[np.inf,0], opt=m.options)

os.system("touch /home/samuroff/hoopoe_paper/toy_model_data/simple_1parfit_e_maxsamples_g1_0.05_g2_0.2_noise8.txt")



gvec=np.linspace(-0.06,0.2,200) 

for i in xrange(30000):
    lvec_noisy0=[]
    noise=np.random.normal(size=32*32).reshape((32,32))*8
    for g1 in gvec:
        im=make_sersic(g1,0.4)
        chi2=((im-g.array+noise)*(im-g.array+noise) ).sum() ; lvec_noisy0.append(-0.5*chi2) ; print g1
    lvec_noisy0=np.array(lvec_noisy0)
    noisy_max = gvec[(lvec_noisy0==lvec_noisy0.max())]

    with open('/home/samuroff/hoopoe_paper/toy_model_data/simple_1parfit_e_maxsamples_g1_0.05_g2_0.2_noise8.txt', "a") as f:
    	f.write('%2.5f \n'%(noisy_max))
    	f.close()

