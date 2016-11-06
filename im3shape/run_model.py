import numpy as np
import galsim
from py3shape import structs
import math
import pylab as plt
import fitsio as fi
import glob, argparse, os
import tools.shapes as s
import py3shape as p3s
from py3shape import utils
from tools.im3shape import basic as i3s
from tools.im3shape import model as toy_model

description = 'Noise free toy model'
parser = argparse.ArgumentParser(description=description, add_help=False)
parser.add_argument('--nrealisations', type=int, action='store', default=1000)
parser.add_argument('--nbins', type=int, action='store', default=12)
args = parser.parse_args()

# Mean neighbour distance (pixels) : 39.79

bins = np.logspace(1,2.8,args.nbins+1)
upper = bins[1:]
lower = bins[:-1]

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

analysis = toy_model.toy_model2(komplexitet=2, random_seed=rank*size)

analysis.setup()

data = s.shapecat()
data.res = fi.FITS("16tiles_results.fits")[1].read()
data.truth = fi.FITS("16tiles_results.fits")[2].read()
ncat = s.shapecat()
ncat.res = fi.FITS("16tiles_neighbours.fits")[1].read()
ncat.truth = fi.FITS("16tiles_neighbours.fits")[2].read()

x=[]
m1=[]
m2=[]
l=[]
u=[]
        
# Cycle over the central flux, as a proxy for SNR
for i, limits in enumerate(zip(lower, upper)):
	if i%size!=rank:
		continue
	snr_min = limits[0]
	snr_max = limits[1]
	print "Level 1 iteration: %d SNR = %3.3f - %3.3f"%(i+1,snr_min, snr_max)

	# Since we're parallelising this, we'll need to deal with SNR bin separately
	# and write out a since file for each
	snr, m11, m22 = analysis.mpi_run(snr_min, snr_max, data, ncat, args.nrealisations)

	x.append(snr)
	l.append(snr_min)
	u.append(snr_max)
	m1.append(m11)
	m2.append(m22)

print "Done all loops"

out=np.vstack((x,l,u,m1,m2))
filename="/home/samuroff/shear_pipeline/toymodel_data-%d.txt"%rank
print "Writing to %s"%filename
np.savetxt(filename, out.T)