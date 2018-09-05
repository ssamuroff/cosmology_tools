import numpy as np 
import itertools
import pylab as plt
plt.style.use('y1a1')

kernel_types = {"e":0, "n":1}
statistics = {"xi+":0, "xi-":1}

columns={True:np.dtype([("bin1", ">i8"), ("bin2", ">i8"), ("corr1", ">i8"), ("corr2", ">i8"), ("q", ">i8"), ("theta1", float), ("bin3", ">i8"), ("bin4", ">i8"), ("corr3", ">i8"), ("corr4", ">i8"), ("r", ">i8"), ("theta2", float), ("cov_gaussian", float), ("cov_nongaussian", float)]),
         False:np.dtype([("bin1", ">i8"), ("bin2", ">i8"), ("corr1", ">i8"), ("corr2", ">i8"), ("ell1", float), ("bin3", ">i8"), ("bin4", ">i8"), ("corr3", ">i8"), ("corr4", ">i8"), ("ell2", float), ("cov_gaussian", float), ("cov_nongaussian", float)])}


class fitscov:
	def __init__(self, filename):
		import fitsio as fi
		print('Reading covariance matrix from %s'%filename)
		self.info = fi.FITS(filename)['covmat'].read_header()
		self.block = fi.FITS(filename)['covmat'].read()
		self.ndim = self.block.shape[0]

		print('Found %dx%d matrix'%(self.ndim,self.ndim))

	def get_eigenvalues(self, force=False):

		if force or (not hasattr(self, 'eigvals')):
			self.eigvals = np.linalg.eigvals(self.block)
		else:
			print('Eigenvalues have already been calculated')

	def get_subblock(self, name):
		order = {'xip':0, 'xim':1, 'gammat':2, 'wtheta':3}

		if isinstance(name, str):
			name1 = name
			name2 = name

		else:
			name1 = name[0]
			name2 = name[1]

		i0 = self.info['STRT_%d'%order[name1]]
		if name1=='wtheta': 
			i1 = -1
		else:
			i1 = self.info['STRT_%d'%(order[name1]+1)]

		j0 = self.info['STRT_%d'%order[name2]]
		if name2=='wtheta': 
			j1 = -1
		else:
			j1 = self.info['STRT_%d'%(order[name2]+1)]

		C = self.block[i0:i1,j0:j1]
		return C

	def diagnose(self, name=None):
		"""Print out some diagnostic statistics either for the whole covariance matrix
		   or for a subblock (defined by a pair of correlation functions)"""

		if name is None:
			C = self.block
		else:
			C = self.get_subblock(name)
			print('Selecting %dx%d subblock'%(C[0].size, C.T[0].size), name)

		symmetric = (C==C.T).all()
		print('Covariance matrix is symmetric:', symmetric)

		try:
			invC = np.linalg.inv(C)
			invertible = True
		except:
			invertible = False
		print('Covariance matrix is invertible:', invertible)

		if name is None:
			self.get_eigenvalues()
			eigvals = self.eigvals
		else:
			eigvals = np.linalg.eigvals(C)
		psd = (eigvals>0).all()
		print('Covariance matrix is positive semidefinite:', psd)
		if not psd:
			nneg = eigvals[(eigvals<0)].size
			print('%d negative eigenvalue(s) : '%nneg, eigvals[(eigvals<0)] )

		det = np.linalg.slogdet(C)
		print('Log determinant: %d'%det[1])

		R = abs(eigvals).max()*1./abs(eigvals).min()
		print('Ratio of largest to smallest eigenvalues: %e/%e = %e'%(abs(eigvals).max(), abs(eigvals).min(), R))

	def show_block(self, name):
		"""Plot out a diagonal subblock of the covariance matrix.
		   Obviously this is obsolete for single-probe covariance matrices."""

		fig = plt.figure(0)
		ax = plt.subplot(111)

		C = self.get_subblock(name)
		
		img = ax.imshow(np.log10(C), origin='lower', interpolation='none')
		plt.show()

	def show(self, partition=True):
		fig = plt.figure(0)
		ax = plt.subplot(111)
		img = ax.imshow(np.log10(self.block), origin='lower', interpolation='none')
		ax.set_aspect('auto')
		plt.xticks([],visible=False)
		plt.yticks([], visible=False)

		if partition:
			labelsy = [r'$\xi_+$', r'$\xi_- $', r'$\gamma_t$', r'$w$']

			i0=0
			i1=0

			for i in range(4):
				plt.axhline(self.info['STRT_%d'%i], color='k')
				plt.axvline(self.info['STRT_%d'%i], color='k')
				i0=i1

				try:
					i1 = self.info['STRT_%d'%(i+1)]
				except:
					i1 = self.block[0].size

				di = (i1-i0)/2
				trans = ax.get_xaxis_transform()
				ax.annotate(labelsy[i], (i0+di-22,1.02), xycoords=trans, fontsize=14)

		plt.show()


	def plot_correlation_matrix(self, savename='/Users/hattifattener/Desktop/multiprobe-multicolour-corrmat-v2.pdf', fontsize=14, cmap='RdBu_r'):

		D = np.sqrt(np.diag(self.block))
		R = np.array([D]*D.size)

		corr = self.block / R / R.T

		fig = plt.figure(0)
		ax = plt.subplot(111)
		img = ax.imshow(corr, origin='lower', interpolation='none', cmap=cmap)
		img.set_clim(-1,1)
		fig.colorbar(img, ax=ax)
		ax.set_aspect('auto')
		plt.xticks([],visible=False)
		plt.yticks([], visible=False)

		labels = [r'$\xi_{+, c1c2}^{ z_{si}, z_{sj}} (\theta)$', r'$\xi_{-, c1c2}^{z_{si},z_{sj}}(\theta)$', r'$\gamma_{t,c2}^{z_{li}, z_{sj}} (\theta)$', r'$w^{z_{li}, z_{lj}} (\theta)$']
		labelsy = [r'$\xi_+$', r'$\xi_- $', r'$\gamma_t$', r'$w$']

		i0=0
		i1=0

		for i in range(4):
			plt.axhline(self.info['STRT_%d'%i], color='k')
			plt.axvline(self.info['STRT_%d'%i], color='k')
			i0=i1
			try:
				i1 = self.info['STRT_%d'%(i+1)]
			except:
				i1 = corr[0].size

			di = (i1-i0)/2
			trans = ax.get_xaxis_transform()
			ax.annotate(labelsy[i], (i0+di-22,1.02), xycoords=trans, fontsize=fontsize)
			ax.annotate(labels[i], (-400, (float(i0+di)/corr[0].size)-0.005), xycoords=trans, fontsize=fontsize)

		plt.savefig(savename)
		plt.close('all')

class covariance_wrapper:
	def __init__(self, covtype, filename, nell=25, nzbins=4, correlations=["ee","ne","nn"], file_type="list", realspace=False):
		if covtype.lower() not in ["gaussian", "hm"]:
			raise ValueError("Please specify covariance type (Gaussian or hm)")

		self.type = covtype.lower()
		self.correlations = correlations

		if file_type=="matrix":
			self.cov = np.loadtxt(filename)

			self.nelem = {"ee":nzbins*nzbins*nell, "ne":nzbins*(nzbins+1)*nell/2, "nn":nzbins*(nzbins+1)*nell/2}

		elif file_type=="list":
			dt = columns[realspace]
			self.cov=np.loadtxt(filename, dtype=dt)
			if realspace:
				xname = "theta"
			else:
				xname="ell"
			setattr(self,xname,np.unique(self.cov["%s1"%xname]))

		self.realspace=realspace

	def extract_all(self, correlations=["ee","ne","nn"], shear_bins=[1,2,3,4], density_bins=[1,2,3, 4], statistics=["cl"], gaussian=True, nongaussian=True, verbosity=1):
		bins = {"e":shear_bins, "n":density_bins}
		covmat={}
		for obs_type1 in statistics:
			for obs_type2 in statistics:
				covmat[(obs_type1,obs_type2)]={}
				if verbosity>0:
					print( obs_type1, obs_type2)
				for c1,c2 in correlations:
					for c3,c4 in correlations:
						covmat[(obs_type1,obs_type2)][(c1,c2,c3,c4)] = {}
						if verbosity>0:
							print( c1,c2, c3, c4)
						bins1 = bins[c1]
						bins2 = bins[c2]
						bins3 = bins[c3]
						bins4 = bins[c4]
				
						for b1 in bins1:
							for b2 in bins2:
								for b3 in bins3:
									for b4 in bins4:
										if verbosity>1:
											print( b1,b2,b3,b4)
										covmat[(obs_type1,obs_type2)][(c1,c2,c3,c4)][(b1,b2,b3,b4)] = self.extract([c1,c2,c3,c4], b1, b2, b3, b4, gaussian=gaussian, nongaussian=nongaussian, realspace=self.realspace, statistic1=obs_type1, statistic2=obs_type2)
										zerolen = covmat[(obs_type1,obs_type2)][(c1,c2,c3,c4)][(b1,b2,b3,b4)].size==0
										if not zerolen:
											continue

										if zerolen:
											if verbosity>2:
												print( "Could not find valid matrix. Will see if we can find it by symmetry arguments.")

											#First try the transpose element
											tmp, zerolen = self.extract([c3,c4,c1,c2], b3, b4, b1, b2, gaussian=gaussian, nongaussian=nongaussian, raise_errors=True, realspace=self.realspace, statistic1=obs_type2, statistic2=obs_type1)
											covmat[(obs_type1,obs_type2)][(c1,c2,c3,c4)][(b1,b2,b3,b4)]=tmp
											if verbosity>2:
												print( "trying transpose block for", b1,b2,b3,b4, obs_type1,obs_type2)
											if tmp.size>0:
												break

											# If the four observable types are all the same all possible permutations of the redshift bins should be identical
#											if (c1==c2==c3==c4) and (obs_type1==obs_type2):
#												for trial_combination in itertools.permutations([b1,b2,b3,b4]):
#													t1,t2,t3,t4 = trial_combination
#													tmp, zerolen = self.extract([c1,c2,c3,c4], t1, t2, t3, t4, gaussian=gaussian, nongaussian=nongaussian, raise_errors=True, realspace=self.realspace, statistic1=obs_type1, statistic2=obs_type2)
#													covmat[(obs_type1,obs_type2)][(c1,c2,c3,c4)][(b1,b2,b3,b4)]=tmp
#													if verbosity>2:
#														print "trying %d %d %d %d for "%(t1,t2,t3,t4), b1,b2,b3,b4, obs_type1,obs_type2
#													if tmp.size>0:
#														break
#												if tmp.size>0:
#														continue

											if (obs_type1!=obs_type2) and zerolen:
												tmp,zerolen = self.extract([c3,c4,c1,c2], b3, b4, b1, b2, gaussian=gaussian, nongaussian=nongaussian, raise_errors=True, realspace=self.realspace, statistic1=obs_type2, statistic2=obs_type1)
												covmat[(obs_type1,obs_type2)][(c1,c2,c3,c4)][(b1,b2,b3,b4)]=tmp
												if verbosity>2:
													print( "reversing datavector order")
												if tmp.size>0:
													break

										for trial_combination in itertools.permutations([(c1,c2,obs_type1,b1,b2),(c3,c4,obs_type2,b3,b4)]):
											(tc1,tc2,tobs_type1,tb1,tb2),(tc3,tc4,tobs_type2,tb3,tb4) = trial_combination

											if "%s%s"%(tc1,tc2) == "%s%s"%(tc3,tc4) and zerolen and (tobs_type1==tobs_type2):
												tmp, zerolen = self.extract([tc1,tc2,tc3,tc4], tb3, tb4, tb1, tb2, gaussian=gaussian, nongaussian=nongaussian, raise_errors=True, realspace=self.realspace, statistic1=tobs_type1, statistic2=tobs_type2)
												covmat[(obs_type1,obs_type2)][(c1,c2,c3,c4)][(b1,b2,b3,b4)]=tmp ; print( "1trying ", tb3, tb4, tb1, tb2)
												if tmp.size>0:
													break

											elif "%s%s"%(tc1,tc2) == "%s%s"%(tc3,tc4) and zerolen and (tobs_type1!=tobs_type2):
												tmp, zerolen = self.extract([tc1,tc2,tc3,tc4], tb3, tb4, tb1, tb2, gaussian=gaussian, nongaussian=nongaussian, raise_errors=True, realspace=self.realspace, statistic1=tobs_type2, statistic2=tobs_type1)
												covmat[(obs_type1,obs_type2)][(c1,c2,c3,c4)][(b1,b2,b3,b4)]=tmp ; print( "2trying ", tb3, tb4, tb1, tb2)
												if tmp.size>0:
													break

											if tc4==tc3 and zerolen:
												tmp, zerolen = self.extract([tc1,tc2,tc3,tc4], tb1, tb2, tb4, tb3, gaussian=gaussian, nongaussian=nongaussian, raise_errors=True, realspace=self.realspace, statistic1=tobs_type1, statistic2=tobs_type2)
												covmat[(obs_type1,obs_type2)][(c1,c2,c3,c4)][(b1,b2,b3,b4)]=tmp ; print( "3trying ", tb1, tb2, tb4, tb3)
												if tmp.size>0:
													break

											elif tc4!=tc3 and zerolen:
												tmp, zerolen = self.extract([tc1,tc2,tc4,tc3], tb1, tb2, tb4, tb3, gaussian=gaussian, nongaussian=nongaussian, raise_errors=True, realspace=self.realspace, statistic1=tobs_type1, statistic2=tobs_type2)
												covmat[(obs_type1,obs_type2)][(c1,c2,c3,c4)][(b1,b2,b3,b4)]=tmp ; print( "4trying ", tb1, tb2, tb4, tb3)
												if tmp.size>0:
													break

											if tc1==tc2 and zerolen:
												tmp, zerolen = self.extract([tc1,tc2,tc3,tc4], tb2, tb1, tb3, tb4, gaussian=gaussian, nongaussian=nongaussian, raise_errors=True, realspace=self.realspace, statistic1=tobs_type1, statistic2=tobs_type2)
												covmat[(obs_type1,obs_type2)][(c1,c2,c3,c4)][(b1,b2,b3,b4)]=tmp ; print( "5trying ", tb2, tb1, tb3, tb4)
												if tmp.size>0:
													break

											elif tc1!=tc2 and zerolen:
												tmp, zerolen = self.extract([tc2,tc1,tc3,tc4], tb2, tb1, tb3, tb4, gaussian=gaussian, nongaussian=nongaussian, raise_errors=True, realspace=self.realspace, statistic1=tobs_type1, statistic2=tobs_type2)
												covmat[(obs_type1,obs_type2)][(c1,c2,c3,c4)][(b1,b2,b3,b4)]=tmp ; print( "6trying ", tb2, tb1, tb3, tb4)
												if tmp.size>0:
													break

											if (tc1==tc2) and (tc4==tc3) and zerolen:
												tmp, zerolen = self.extract([tc1,tc2,tc3,tc4], tb2, tb1, tb4, tb3, gaussian=gaussian, nongaussian=nongaussian, raise_errors=True, realspace=self.realspace, statistic1=tobs_type1, statistic2=tobs_type2)
												covmat[(obs_type1,obs_type2)][(c1,c2,c3,c4)][(b1,b2,b3,b4)]=tmp ; print( "7trying ", tb2, tb1, tb4, tb3)
												if tmp.size>0: 
													break

										if tmp.size>0:
											continue
											
        									
        
										if zerolen:
											print( "No covariance block could be found.")
											import pdb ; pdb.set_trace()
		return covmat


	def extract(self, corr, bin1, bin2, bin3, bin4, verbose=False, gaussian=True, nongaussian=True, raise_errors=False, realspace=False, statistic1="xi+", statistic2="xi+"):
		k1 = kernel_types[corr[0]]
		k2 = kernel_types[corr[1]]
		k3 = kernel_types[corr[2]]
		k4 = kernel_types[corr[3]]

		if verbose:
			print( "kernel types: ",k1,k2,k3,k4)

		sel = (self.cov["corr1"]==k1) & (self.cov["corr2"]==k2) & (self.cov["corr3"]==k3) & (self.cov["corr4"]==k4) & (self.cov["bin1"]==bin1) & (self.cov["bin2"]==bin2) & (self.cov["bin3"]==bin3) & (self.cov["bin4"]==bin4)

		cov = self.cov[sel]

		if realspace:
			x = np.unique(cov["theta1"])
			xname="theta"
		else:
			x = np.unique(cov["ell1"])
			xname="ell"
		nx = len(x)

		m = np.zeros((nx,nx))

		# Choose the relevant observable
		if realspace:
			sel3 = ((cov["q"]==statistics[statistic1]) & (cov["r"]==statistics[statistic2]))
			reverse_bins=False
			if len(cov["q"][sel3])==0:
				sel3 = ((cov["q"]==statistics[statistic2]) & (cov["r"]==statistics[statistic1]))
				reverse_bins=True
		else:
			sel3 = np.ones_like(cov).astype(bool)



		for i1,x1 in enumerate(x):
			for i2 in xrange(i1,nx):
				x2 = x[i2]
				#if reverse_bins:
					#sel2 = (cov["%s1"%xname]==x2) & (cov["%s2"%xname]==x1)
				#else:
				sel2 = (cov["%s1"%xname]==x1) & (cov["%s2"%xname]==x2)
				
				# Pull out Gaussian and non Gaussian contributions
				if gaussian:
					cg = cov["cov_gaussian"][sel2 & sel3]
				else:
					cg=0
				if nongaussian:
					cng = cov["cov_nongaussian"][sel2 & sel3]
				else:
					cng=0
				try:
					m[i2,i1] = cg+cng
				except:
					import pdb; pdb.set_trace()
				m[i1,i2] = cg+cng
				if verbose:
					print( i1,i2)


		zero_length = len(m)==0
		if not raise_errors:
			return m
		else:
			return m, zero_length

	def plot(self, log=True, submatrix=None):
		import pylab as plt

		if submatrix:
			c = getattr(self,submatrix)
		else:
			c = self.cov
		if log:
			c = np.log(c)

		plt.matshow(c)
		plt.colorbar()
		plt.title(self.type)

	def slice(self):
		i=0
		j = 0
		for corr1 in self.correlations:
			for corr2 in self.correlations:
				print( "Extracting sub matrix: %s%s"%(corr1,corr2))
				nm1 = self.nelem[corr1]
				nm2 = self.nelem[corr2]
				sub = self.cov[j:j+nm1,i:i+nm2]

				setattr(self, "%s%s"%(corr1,corr2), sub)
				i+=nm2

			j+=nm1
			i=0

	def format(self):
		ctypes=["ee","ne","nn"]
		newmat = np.zeros_like(self.cov)
		i=0
		j=0

		for corr1 in ctypes:
			for corr2 in ctypes:
				nm1 = self.nelem[corr1]
				nm2 = self.nelem[corr2]
				print( "Inserting %s%s (%d-%d,%d-%d)"%(corr1,corr2,j,j+nm1,i,i+nm2))
				c = getattr(self,"%s%s"%(corr1,corr2))

				newmat[j:j+nm1,i:i+nm2] = c

				i+=nm2

			j+=nm1
			i=0

	def write(self,filename):
		np.savetxt(filename,self.corr)


order={'xip':(0,1), 'xim':(1,2), 'gammat':(2,3), 'wtheta':(3,-1)}
def compute_snr(fits, corr='all'):

	if (corr=='all'):
		# Assemble the 3x2pt datavector
		dvec = np.concatenate([fits[c]['VALUE'].read() for c in ['xip','xim','gammat','wtheta']])

		# And the covaraiance matrix
		C = fits['COVMAT'].read()
	else:
		dvec = fits[corr]['VALUE'].read()

		C = fits['COVMAT'].read()
		hdr = fits['COVMAT'].read_header()
		i0 = hdr['STRT_%d'%order[corr][0]]
		if (order[corr][1]==-1):
			C = C[i0:,i0:]
		else:
			i1 = hdr['STRT_%d'%order[corr][1]]
			C = C[i0:i1,i0:i1]
		

	print('Datavector contains %d elements.'%len(dvec))
	print('Inverting covariance matrix.')
	Cinv = np.linalg.inv(C)

	snr = np.sqrt(np.dot(dvec, np.dot(Cinv,dvec) ))
	print('S/R: %3.3f'%snr)

	return snr







