import numpy as np
import fitsio as fi
import argparse
import yaml
import scipy.special as sps


description = 'Covariance building code.'
parser = argparse.ArgumentParser(description=description, add_help=True)
parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3), help='integer verbosity level: min=0, max=3 [default=2]')
parser.add_argument('-o', '--order', type=str, action='store', help='FITS file from which to obtain data order.')
parser.add_argument('-t', '--theory', type=str, action='store', help='Directory containing angular power predictions from theory.')
parser.add_argument('-c', '--config', type=str, action='store', help='Configuration file containing survey specific parameter values.')
args = parser.parse_args()

datavec = fi.FITS(args.order)

params = yaml.load(open(args.config))

area = (np.pi * params['area']/180.)
fsky = area/4./np.pi

bessel_kernels = {'xip':0, 'xim':4, 'gammat':2}
corr_order = {'xip':0, 'xim':1, 'gammat':2, 'wtheta':3}

class matrix:
	def __init__(self, datavec, theory):
		print 'Initialised covariance matrix'
		self.base = theory
		# Bad, bad, bad. I know. This is bad.
		# Will fix this at some point
		naxis = 2080 + 700 + 100
		self.block = np.zeros((2100,2100)) + 1e-30
	def choose_spectra(self, i, j, k, l, correlation1, correlation2):
		# For the (xi xi) part all of the required angular spectra are GG
		# The red-blue information is subsumed into the index 
		if ('xi' in correlation1) and ('xi' in correlation2):
			path = '%s/multicolour-3x2pt/shear_cl/'%(self.base) + 'bin_%d_%d.txt'
			C1 = np.genfromtxt(path%(k,j))
			C2 = np.genfromtxt(path%(l,i))
			C3 = np.genfromtxt(path%(k,i))
			C4 = np.genfromtxt(path%(l,j))
		# For the (gammat gammat) all of the GG, gg and gG Cls are needed 
		elif (correlation1==correlation2) and ('gammat' in correlation1):
			path_ss = '%s/multicolour-3x2pt/shear_cl/'%(self.base) + 'bin_%d_%d.txt'
			path_gs = '%s/multicolour-3x2pt/galaxy_shear_cl/'%(self.base) + 'bin_%d_%d.txt'
			path_gg = '%s/multicolour-3x2pt/galaxy_cl/'%(self.base) + 'bin_%d_%d.txt'
			C1 = np.genfromtxt(path_gs%(k,j))
			C2 = np.genfromtxt(path_gs%(l,i))
			C3 = np.genfromtxt(path_ss%(k,i))
			C4 = np.genfromtxt(path_gg%(l,j))
			# Shot noise term
			if (j==l):
				C4 += 1/params['lens_n_gal'][j]
		else:
			path_ss = '%s/multicolour-3x2pt/shear_cl/'%(self.base) + 'bin_%d_%d.txt'
			path_gs = '%s/multicolour-3x2pt/galaxy_shear_cl/'%(self.base) + 'bin_%d_%d.txt'
			C1 = np.genfromtxt(path_gs%(k,j))
			C2 = np.genfromtxt(path_ss%(l,i))
			C3 = np.genfromtxt(path_ss%(k,i))
			C4 = np.genfromtxt(path_ss%(l,j))

		x = np.genfromtxt('%s/multicolour-3x2pt/shear_cl/ell.txt'%self.base)

		return x, C1, C2, C3, C4
	def find_element(self, t1, t2, i, j, k, l, gtype, correlation1, correlation2):
		if not hasattr(self, gtype):
			path = args.order.replace('multicolour', gtype)
			setattr(self, gtype, fi.FITS(path)['covmat'.upper()].read())
			setattr(self, '%s_info'%gtype, fi.FITS(path)['covmat'.upper()].read_header())
			setattr(self, '%s_%s'%(gtype, 'xip'), fi.FITS(path)['xip'].read())
			setattr(self, '%s_%s'%(gtype, 'xim'), fi.FITS(path)['xim'].read())
			setattr(self, '%s_%s'%(gtype, 'gammat'), fi.FITS(path)['gammat'].read())
			setattr(self, '%s_%s'%(gtype, 'wtheta'), fi.FITS(path)['wtheta'].read())

			
		info = getattr(self, '%s_info'%gtype)
		cov = getattr(self, gtype)
		ord1 = getattr(self, '%s_%s'%(gtype, correlation1))
		ord2 = getattr(self, '%s_%s'%(gtype, correlation2))

		# Look up the correct position along each axis
		index1 = np.argwhere((ord1['BIN1']==i) & (ord1['BIN2']==j) & (ord1['ANGBIN']==t1))[0,0]
		index2 = np.argwhere((ord2['BIN1']==k) & (ord2['BIN2']==l) & (ord2['ANGBIN']==t2))[0,0]

		return cov[index2, index1]

	def export(self):
		print 'Saving covariance matrix'
		outfits = fi.FITS('covmat.fits', 'rw')
		outfits.write(self.block)
		outfits[-1].write_key('EXTNAME', 'COVMAT')
		outfits.close()

	def interpret_bin_indices(self, i, j, corr):
		if (corr=='wtheta'):
			return 'lens', 'lens',i, j
		elif (corr=='gammat'):
			if j>4:
				j-=4
				source_type = 'early'
			else:
				source_type = 'late'
			return 'lens', source_type, i, j

		elif ('xi' in corr):
			if j>4:
				j-=4
				source_type2 = 'early'
			else:
				source_type2 = 'late'
			if i>4:
				i-=4
				source_type1 = 'early'
			else:
				source_type1 = 'late'
			return source_type1, source_type2, i, j


cov = matrix(datavec, args.theory)
a,b = 0,0

for correlation1 in ['xip', 'xim', 'gammat', 'wtheta']:
	for (itheta1, theta1, i,j) in zip(datavec[correlation1]['ANGBIN'][:], datavec[correlation1]['ANG'][:], datavec[correlation1]['BIN1'][:], datavec[correlation1]['BIN2'][:]):

			row = []
			for correlation2 in ['xip', 'xim', 'gammat', 'wtheta']:
				for (itheta2, theta2, k,l) in zip(datavec[correlation2]['ANGBIN'][:], datavec[correlation2]['ANG'][:], datavec[correlation2]['BIN1'][:], datavec[correlation2]['BIN2'][:]):

						g1, g2, i0, j0 = cov.interpret_bin_indices(i, j, correlation1)
						g3, g4, k0, l0 = cov.interpret_bin_indices(k, l, correlation2)

						if args.verbosity>1:
							print '%s %s'%(correlation1, correlation2)
							print 'Colour bins: (%s,%s),(%s,%s)'%(g1,g2,g3,g4)
							print 'Redshift bins: (%d,%d),(%d,%d)'%(i0,j0,k0,l0)
							print 'Angular bins: (%d,%d)'%(itheta1,itheta2)

						source_types = [g for g in [g1,g2,g3,g4] if g!='lens']

						if (len(np.unique(source_types))>1):
							try:
								ell, C1, C2, C3, C4 = cov.choose_spectra(i, j, k, l, correlation1, correlation2)
							except:
								cov.block[b,a] = cov.block[a,b]
								continue

								

							n1 = bessel_kernels[correlation1]
							n2 = bessel_kernels[correlation2]
							J1 = sps.jn(n1, ell*theta1)
							J2 = sps.jn(n2, ell*theta2)

							K = ell * J1 * J2 * (C1*C2 + C3*C4)

							val = np.trapz(K, ell) / (8*np.pi*np.pi) / fsky
						else:
							val = cov.find_element(itheta1, itheta2, i, j, k, l, source_types[0], correlation1, correlation2)

						cov.block[b,a] = val
						a+=1
			b+=1
			a=0



import pdb ; pdb.set_trace()
cov.export()






				