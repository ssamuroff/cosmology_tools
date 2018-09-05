# See www.slac.stanford.edu/~jderose/buzzard/buzzard_v1.1.html#y1-bright

import glob, treecorr
import scipy as sp
import numpy as np 
import pyfits

class buzzard:
	cosmology = {"Omega_m" : 0.286,
		"Omega_de" : 0.714,
		"ns" : 0.96,
		"h" : 0.7}

	def __init__(self, nfiles, mode='obs'):
		if mode=='obs': self.files = glob.glob("Buzzard*.fit")
		self.nfil = nfiles
		self.cat = pyfits.getdata(self.files[0])
		for i in range(1,nfiles): self.cat = np.concatenate((self.cat,pyfits.getdata(self.files[i]))) 
		self.ngal = len(self.cat)
		print 'Got %2.2f M galaxies.' %(self.ngal/1.0e6)

	def get_sky_patch(self):
		ra_max = self.cat['RA'].max()
		ra_min = self.cat['RA'].min()
		ra0 = (ra_max + ra_min) /2.
		dra = 1.0
		sel = (self.cat['RA'] > (ra0-dra) ) & (self.cat['RA'] < (ra0+dra) )
		if len(self.cat[sel])==0:
			sel = (self.cat['RA'] < ra0 )
			self.cat = self.cat[sel]
		print 'Selected sky area contains %f M galaxies.'%(len(self.cat)/1.e6)

	@staticmethod
	def info():
		print 'Please refer to '
		print 'www.slac.stanford.edu/~jderose/buzzard/buzzard_v1.1.html#y1-bright'

		print cosmology
		

	@staticmethod
	def get_zbins(cat, nzbin):
		edges = [0]
		for i in range(1,nzbin+1):
	     		edges += [sp.percentile(cat['PHOTOZ_GAUSSIAN'], i*100./nzbin)]
	     		print 'bin %d'%i

		return edges

	@staticmethod
	def save_correlation(bz, ctype, nz, outpath):
		edges =  bz.get_zbins(bz.cat, nz)

		out = {'z_bin_edges': edges, 'nzbins': nz}

		for i in range (1,nz+1):
			if ctype[0] == ctype[1]: jlim = i
			else: jlim = 0
			for j in range(jlim,nz+1):
				exec 'corr =  bz.%s%d%d' %(ctype,i,j)
				if ctype == 'GG': 
					corr1 = 'xip' ; corr2 = 'xim'
					var = 'varxi' ; npair = 'npairs' ; w = 'weight'
				exec "out['bin_%d_%d'] = ( corr.%s , corr.%s, corr.%s, corr.%s, corr.%s)" %(i, j, corr1, corr2, var, npair, w)
				if (i,j) == (1,1): out['theta'] = np.exp(corr.logr) 
		

	def setup_correlations(self, corrtypes=['GG'], nzbin=6, ntbin=12, thetalim=(0.1,40.)):
		self.nzbin = nzbin
		self.ntbin = ntbin

		edges = self.get_zbins(self.cat, nzbin)
		s="cat%d=treecorr.Catalog(ra=self.cat[sel]['RA'], dec=self.cat[sel]['DEC'], ra_units='degrees', dec_units='degrees', g1=self.cat[sel]['EPSILON1'], g2=self.cat[sel]['EPSILON2'])"
		s2="self.%s%d%d=treecorr.%sCorrelation(min_sep=%f, max_sep= %f, sep_units='arcmin', bin_size=%f/(%f-%f) )"
		for i in range (1,nzbin+1):
			sel = (self.cat['PHOTOZ_GAUSSIAN']>edges[i-1]) & (self.cat['PHOTOZ_GAUSSIAN']<edges[i])
			exec s%i

		for i in range (1,nzbin+1):
			jlim = i
			if 'NG' in corrtypes: jlim = 0
	     		for j in range(jlim,nzbin+1):
				for corr in corrtypes:
					exec s2%(corr, i,j, corr, thetalim[0], thetalim[1], ntbin, thetalim[1], thetalim[0])
					exec "self.%s%d%d.process(cat%d, cat%d)"%(corr,i,j,i,j)
					print corr, i, j

	def get_n_of_z(self,zlim=(0.,2.5), nzbin=6):
		dat2=pyfits.getdata(self.files[0])

		edges = self.get_zbins()
		
		s3="self.nz_%d+=np.exp(-1.*(zf-zpoint[j])**2 /2/(sigp[j]*sigp[j]))"
		z = dat2['PHOTOZ_GAUSSIAN']
		zf = np.linspace(zlim[0],zlim[1],400)
		self.pz= {'z': zf, 'z_bin_edges': edges}
		sig = 0.03/(1.+z)
		for i in range(1,nzbin+1):
			self.pz['bin_%d'%i]=np.zeros_like(zf)	     
			sel = (z>edges[i-1]) & (z<edges[i])
			sigp = sig[sel]
			zpoint = z[sel]
			for j in range (len(zpoint)):
				self.pz['bin_%d'%i]+=np.exp(-1.*(zf-zpoint[j])**2 /2/(sigp[j]*sigp[j]))
			print i,j

	def get_luminosity_function(self, band='r', nbin = 50):
		ind = {'g': 0, 'r': 1, 'i': 2, 'z': 3}
		r = self.cat['AMAG'].T[ind[band]]

		redshift = self.cat['Z']

		arg = -1.*r*2./5.
		arg += np.log10(4.*np.pi)*2./5.
		arg += 5.
		L = 10**arg

		Lmin = L.min()
		Lmax = L.max()

		zmin = 0.
		zmax=np.ceil(redshift.max())

		lbins = np.logspace(np.log10(Lmin), np.log10(Lmax), nbin)
		zbins = np.linspace(zmin, zmax, nbin)

		phi = np.zeros((nbin,nbin))
		l = np.zeros((nbin,nbin))
		V = np.zeros((nbin,nbin))

		dA = 1000.*np.pi*np.pi/180./180.
		c = 3.0e8
		for i in xrange(nbin-1):
			sel = (redshift>zbins[i]) & (redshift<zbins[i+1])
			for j in xrange(nbin-1):
				sel2 = (L[sel]>lbins[j]) & (L[sel]<lbins[j+1])
				N = len(L[sel][sel2])
				dV = 0.
				if N!=0:
					dz = redshift[sel][sel2].max() - redshift[sel][sel2].min()
					z0 = np.mean(redshift[sel][sel2])
					z = np.linspace(redshift[sel][sel2].min(), redshift[sel][sel2].max(), 1000)
					dx = np.trapz(np.sqrt(self.cosmology['Omega_m']*(1+z)**-3 + self.cosmology['Omega_de']), x=z  )
					z = np.linspace(0., redshift[sel][sel2].max(), 1000)
					x0 = np.trapz(np.sqrt(self.cosmology['Omega_m']*(1+z)**-3 + self.cosmology['Omega_de']), x=z )
					# Comoving volume element in this bin
					dV = x0 * x0 * dA * dx
					# The volume is in units of h^3
					dV *= c*c*c/(100*100*100)

				
					phi[i,j] = N / dV 
					l[i,j] = np.mean(L[sel][sel2])
					V[i,j] = dV
				print 'phi[%d,%d] = %e (%d)' %(i, j, phi[i,j], N)

		return phi, l, V, zbins, lbins 

	def red_cuts(self):
		g,r,i,z,Y = self.cat["OMAG"].T

		# Cuts designed to give a red sample for BAO measurements
		# Not the same as redmagic, which includes more  than colour cuts 
		colour_cuts = (i>17.5) & (i<22.) & (g-r > -1) & (g-r < 3) 
		colour_cuts = colour_cuts & (r-i > -1) & (r-i < 2.5)
		colour_cuts = colour_cuts & (i-z> -1) & (i-z < 2.)

		print '%f M galaxies pass all colour cuts.'%(len(self.cat[colour_cuts])/1.e6)

		return colour_cuts
