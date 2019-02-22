import fitsio as fi
import numpy as np
import os

redshifts = np.array([0.062,0.300,0.625,1.000])
names = ['wgg','wgp','wpp']
covfile = 'shape_noise-covmat.txt'

os.system('rm 2pt_ia_blindtheory_sncov.fits')
outfits = fi.FITS('2pt_ia_blindtheory_sncov.fits','rw')

for hduname in names:
	print('Processing %s'%hduname)
	out = {'VALUE':[],'SEP':[],'BIN':[]}
	for i,z in enumerate(redshifts):
		filename = '%s_%3.3f.txt'%(hduname,z)
		r,w = np.loadtxt(filename).T
		index = [i]*len(r)
		out['BIN'].append(index)
		out['VALUE'].append(w)
		out['SEP'].append(r)

	for col in out.keys():
		out[col] = np.concatenate(out[col])

	outfits.write(out)
	outfits[-1].write_key('EXTNAME', hduname)
	outfits[-1].write_key('2PTDATA', True)
	outfits[-1].write_key('N_R', len(r))
	outfits[-1].write_key('TUNIT2', 'h-1 Mpc')

print('Including covariance matrix: %s'%covfile)
cov = np.loadtxt(covfile)
outfits.write(cov)
outfits[-1].write_key('EXTNAME', 'COVMAT')
outfits[-1].write_key('NREAL', 1111)

outfits.close()




