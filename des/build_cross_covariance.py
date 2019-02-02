import numpy as np
import fitsio as fi
import argparse
import yaml
import os
import scipy.special as sps


description = 'Covariance building code.'
parser = argparse.ArgumentParser(description=description, add_help=True)
parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3), help='integer verbosity level: min=0, max=3 [default=2]')
parser.add_argument('-o', '--order', type=str, action='store', help='FITS file from which to obtain data order.')
parser.add_argument('-t', '--theory', type=str, action='store', help='Directory containing angular power predictions from theory.')
parser.add_argument('-c', '--config', type=str, action='store', help='Configuration file containing survey specific parameter values.')
parser.add_argument('--outdir', type=str, default='/home/ssamurof/cosmosis/cov', action='store', help='Directory to which to save the output covariance matrix.')
parser.add_argument('--mpi', action='store_true', help='Parallelize job?')
args = parser.parse_args()

if args.mpi:
	import mpi4py.MPI
	rank = mpi4py.MPI.COMM_WORLD.Get_rank()
	size = mpi4py.MPI.COMM_WORLD.Get_size()
	print 'Setting up MPI with %d threads (thread %d).'%(size, rank)
else:
	rank = 0
	size = 1

datavec = fi.FITS(args.order)

params = yaml.load(open(args.config))

area = (np.pi * np.pi * params['area']/180./180.)
fsky = area/4./np.pi

bessel_kernels = {'xip':0, 'xim':4, 'gammat':2, 'wtheta':0}
corr_order = {'xip':0, 'xim':1, 'gammat':2, 'wtheta':3}

def get_cl(path, i,j, ctype='ss', noise=True):

	# Decide whether a noise term is needed 
	k = (60.*180./np.pi)**2

	if (not noise):
		N = 0.0
	elif (ctype=='ss') and (i==j):
		ngal = params['source_n_gal'][j-1] * k
		N = params['sigma_e'][j-1]*params['sigma_e'][j-1]/ngal/2
	elif (ctype=='gg') and (i==j):
		ngal = params['lens_n_gal'][j-1] * k
		N = 1/ngal
	else:
		N = 0.0

	if os.path.exists(path%(i,j)):
		return np.genfromtxt(path%(i,j)) + N
	else:
		return np.genfromtxt(path%(j,i)) + N


class matrix:
	def __init__(self, datavec, theory):
		print 'Initialised covariance matrix'
		self.base = theory
		# Use the FITS file to decide the shape of the covariance matrix
		hdr = datavec['covmat'].read_header()
		naxis = hdr['naxis1']
		# Set up an empty square array to which to write 
		self.block = np.zeros((naxis,naxis)) + 1e-30

	def choose_spectra(self, i, j, k, l, correlation1, correlation2):
		# Now we need to determine the four relevant Cl terms to use 
		# for the cosmic variance part of the covariance 
		# This is the fiddly bit

		# eq : Cov[C_ab^ij, C_cd^kl] \prop (C_ac^ik x C_bd^jl + C_ad^il x C_bc^jk) 

		path_gs = '%s/multicolour-3x2pt/galaxy_shear_cl/'%(self.base) + 'bin_%d_%d.txt'
		path_ss = '%s/multicolour-3x2pt/shear_cl/'%(self.base) + 'bin_%d_%d.txt'
		path_gg = '%s/multicolour-3x2pt/galaxy_cl/'%(self.base) + 'bin_%d_%d.txt'

		# For the (xi xi) part all of the required angular spectra are GG
		# But the calculation should not include a shape noise term in the <xip xim> cross correlations 
		if ('xi' in correlation1) and ('xi' in correlation2):
			include_noise = (correlation1==correlation2)
			C1 = get_cl(path_ss, k, j, ctype='ss', noise=include_noise)
			C2 = get_cl(path_ss, l, i, ctype='ss', noise=include_noise)
			C3 = get_cl(path_ss, k, i, ctype='ss', noise=include_noise)
			C4 = get_cl(path_ss, l, j, ctype='ss', noise=include_noise)

		elif  (('gammat' in correlation2) and ('xi' in correlation1)):
			C1 = get_cl(path_gs, k, i, ctype='gs')
			C2 = get_cl(path_ss, j, l, ctype='ss')
			C3 = get_cl(path_ss, i, l, ctype='ss')
			C4 = get_cl(path_gs, k, j, ctype='gs')

		elif  (('wtheta' in correlation2) and ('xi' in correlation1)):
			C1 = get_cl(path_gs, k, i, ctype='gs')
			C2 = get_cl(path_gs, l, j, ctype='gs')
			C3 = get_cl(path_gs, l, i, ctype='gs')
			C4 = get_cl(path_gs, k, j, ctype='gs')

		elif  (('gammat' in correlation2) and ('gammat' in correlation1)):
			C1 = get_cl(path_gg, k, i, ctype='gg')
			C2 = get_cl(path_ss, l, j, ctype='ss')
			C3 = get_cl(path_gs, l, i, ctype='gs')
			C4 = get_cl(path_gs, k, j, ctype='gs')

		elif  (('wtheta' in correlation2) and ('gammat' in correlation1)):
			C1 = get_cl(path_gg, k, i, ctype='gg')
			C2 = get_cl(path_gs, l, j, ctype='gs')
			C3 = get_cl(path_gg, l, i, ctype='gg')
			C4 = get_cl(path_gs, k, j, ctype='gs')

		elif  (('gammat' in correlation2) and ('wtheta' in correlation1)):
			C1 = get_cl(path_gg, k, i, ctype='gg')
			C2 = get_cl(path_gs, j, l, ctype='gs')
			C3 = get_cl(path_gs, i, l, ctype='gs')
			C4 = get_cl(path_gg, k, j, ctype='gg')

		elif  (('wtheta' in correlation2) and ('wtheta' in correlation1)):
			C1 = get_cl(path_gg, k, i, ctype='gg')
			C2 = get_cl(path_gg, l, j, ctype='gg')
			C3 = get_cl(path_gg, l, i, ctype='gg')
			C4 = get_cl(path_gg, k, j, ctype='gg')

		x = np.genfromtxt('%s/multicolour-3x2pt/shear_cl/ell.txt'%self.base)

		return x, C1, C2, C3, C4

	def find_element(self, t1, t2, i, j, k, l, gtype, correlation1, correlation2):
		if not hasattr(self, gtype):
			path = args.order.replace('multicolour', gtype).replace('_data','').replace('-v7', '')
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
		try:
			index1 = np.argwhere((ord1['BIN1']==i) & (ord1['BIN2']==j) & (ord1['ANGBIN']==t1))[0,0]
		except:
			import pdb ; pdb.set_trace()
			if ('xi' in correlation1) or ('wtheta' in correlation1):
				index1 = np.argwhere((ord1['BIN1']==j) & (ord1['BIN2']==i) & (ord1['ANGBIN']==t1))[0,0]
			else:
				import pdb ; pdb.set_trace()

                #if ((g1,g2,g3,g4)==('early','early','lens','early')) and ('xi' in correlation1) and ('gammat' in correlation2) and ((i,j,k,l)==(1,1,3,1)): import pdb ; pdb.set_trace()

		try:
			index2 = np.argwhere((ord2['BIN1']==k) & (ord2['BIN2']==l) & (ord2['ANGBIN']==t2))[0,0]
		except:
			if ('xi' in correlation2) or ('wtheta' in correlation2):
				index2 = np.argwhere((ord2['BIN1']==l) & (ord2['BIN2']==k) & (ord2['ANGBIN']==t2))[0,0]
			else:
				raise

		start1 = info['STRT_%d'%corr_order[correlation1]]
		start2 = info['STRT_%d'%corr_order[correlation2]]

		return cov[start2+index2, start1+index1]

	def export(self, suffix=''):
		print 'Saving covariance matrix'
		outfits = fi.FITS('/home/samuroff/covmat-%d.fits'%suffix, 'rw')
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

outfile = open('%s/cov-out-%s-v0.txt'%(args.outdir, rank), 'wa')
ielem=0

for correlation1 in ['xip','xim','gammat','wtheta']:
	for (itheta1, theta1, i,j) in zip(datavec[correlation1]['ANGBIN'][:], datavec[correlation1]['ANG'][:], datavec[correlation1]['BIN1'][:], datavec[correlation1]['BIN2'][:]):

		for correlation2 in ['xip','xim','gammat','wtheta']:
			for (itheta2, theta2, k, l) in zip(datavec[correlation2]['ANGBIN'][:], datavec[correlation2]['ANG'][:], datavec[correlation2]['BIN1'][:], datavec[correlation2]['BIN2'][:]):

				# MPI handling
				if (ielem%size!=rank):
					ielem+=1
					a+=1
					continue

				# If we can find this element by symmetry then do so
				if (a!=b) and (abs(cov.block[a,b])>1e-30):
					val = 1.0 * cov.block[a,b]
					cov.block[b,a] = val
					outfile.write('%d %d %e \n'%(a, b, val))
					a+=1
					ielem+=1
					continue

				g1, g2, i0, j0 = cov.interpret_bin_indices(i, j, correlation1)
				g3, g4, k0, l0 = cov.interpret_bin_indices(k, l, correlation2) 

				if args.verbosity>1:
					print '%s %s'%(correlation1, correlation2)
					print 'Colour bins: (%s,%s),(%s,%s)'%(g1,g2,g3,g4)
					print 'Redshift bins: (%d,%d),(%d,%d)'%(i0,j0,k0,l0)
					print 'Angular bins: (%d,%d)'%(itheta1,itheta2)

				source_types = [g for g in [g1,g2,g3,g4] if g!='lens']

				if (len(np.unique(source_types))>1):
					#val = 1e-15

					ell, C1, C2, C3, C4 = cov.choose_spectra(i, j, k, l, correlation1, correlation2)
					n1 = bessel_kernels[correlation1]
					n2 = bessel_kernels[correlation2]
					J1 = sps.jn(n1, ell*theta1)
					J2 = sps.jn(n2, ell*theta2)

					K = ell * J1 * J2 * (C1*C2 + C3*C4)

					val = np.trapz(K, ell) / (8*np.pi*np.pi) / fsky
				else:
					if len(source_types)==0:
						gt = 'late'
					else:
						gt = source_types[0]

					val = cov.find_element(itheta1, itheta2, i0, j0, k0, l0, gt, correlation1, correlation2)
					#except:
					#	ell, C1, C2, C3, C4 = cov.choose_spectra(i, j, k, l, correlation1, correlation2)
					#	n1 = bessel_kernels[correlation1]
					#	n2 = bessel_kernels[correlation2]
					#	J1 = sps.jn(n1, ell*theta1)
					#	J2 = sps.jn(n2, ell*theta2)
					#	K = ell * J1 * J2 * (C1*C2 + C3*C4)
					#	val = np.trapz(K, ell) / (8*np.pi*np.pi) / fsky
                                            
					#val = 5e-14

				cov.block[b,a] = val
				print val 
				#if val==0.0: import pdb ; pdb.set_trace()
				
				outfile.write('%d %d %e \n'%(a, b, val))
				a+=1
				ielem+=1
		b+=1
		a=0

outfile.close()
cov.export(suffix=rank)






				
