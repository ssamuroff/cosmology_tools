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
from tools.im3shape import basic as i3s

description = 'Noise free toy model'
parser = argparse.ArgumentParser(description=description, add_help=False)
parser.add_argument('-p', '--parameter', type=str, action='store', default="psf_size")
args = parser.parse_args()
param_name=args.parameter

bins={"neighbour_size":np.array([1,1.5,2.5,3,3.5,4,4.5,5]), "psf_size":np.array([1,1.5,2.5,3,3.5,4,4.5,5]), "psf_ellipticity":np.linspace(-1,1,9), "neighbour_flux":np.linspace(1000,3000,9), "central_flux":np.linspace(1000,3000,9), "dneigh":np.linspace(5,30,9)}

print "Will calculate alpha samples"
print "-------------------------------------------"

q = bins[param_name]
print "parameter:", param_name
print "sampling:", q

model = i3s.toy_model()

x,alpha = model.sample_alpha(param_name, q)

out=np.vstack((x,alpha))
filename="/home/samuroff/shear_pipeline/end-to-end/%s_toymodel_data.txt"%param_name
print "Writing to %s"%filename
np.savetxt(filename, out.T)