import numpy as np
import scipy as sp
import scipy.interpolate as spi
import scipy.spatial as sps
import astropy.io.fits as pyfits

#import matplotlib.colors
#import matplotlib
import math
#matplotlib.rcParams['font.family']='serif'
#matplotlib.rcParams['font.size']=14
#matplotlib.rcParams['legend.fontsize']=14
#matplotlib.rcParams['xtick.major.size'] = 10.0
#matplotlib.rcParams['ytick.major.size'] = 10.0

colours=['midnightblue','forestgreen','pink', 'purple', "lightblue"]
constants={"c" : 299792458.0}

class calculator:
	def load_from_txt(self, filename=None):
		self.pz=[]
		if not filename:
			filename=self.files
		dat = np.loadtxt(filename[0]).T
		self.z = dat[0]
		for p in dat[1:]:
		    self.pz.append(p)

		print "Loaded n(z) in %d tomographic bins from %s"%(len(self.pz), filename)


class n_z_calc:
	def __init__(self, filename, mode=None):
		self.mode = mode
		if (mode=='text'):
			dat = np.loadtxt(filename)

			ra = dat.T[1]
			dec = dat.T[2]
			z = dat.T[3]
			dz = dat.T[4]
			self.data = {'ra': ra , 'dec': dec, 'point_z': z, 'sigma_z': dz}
		if (mode=='skynet'):
			self.hdus = pyfits.open(filename)
			self.data = {}
		else: 
			print 'Unrecognised file format. Please specify mode = txt or skynet.'
			return None

	def select_tiles(self, tilelist, fmt='txt'):
		if fmt=='txt':
			tiles = np.loadtxt(tilelist)
		elif fmt=='npy':
			tiles = np.load(tilelist)
		else:
			print 'Unrecognised file format.'
			return None
		return self.get_mask(self.data['coadd_objects_id'], tiles)
		
	@staticmethod
	def get_mask(coadd_tiles, selection):
		sel = np.in1d(coadd_tiles,selection)
		print 'Using %d of %d galaxies.'%(len(coadd_tiles[sel]), len(coadd_tiles))
		return sel
			
	def extract_columns(self, nsamp, bin_centres=True, ids=True, means=True, modes=True, pz_samples=True, stdv=True):
		if self.mode=='skynet':
			if bin_centres: 
				zc = [zp[0] for zp in self.hdus[1].data]
				self.data.update({'bin_centres':np.array(zc)})
				edges = (np.array(zc)[1:]+np.array(zc)[:-1])/2.
				lim1 = np.array(zc)[0] - (np.array(zc)[1]-np.array(zc)[0])
				if lim1<0.0: lim1=0.0
				lim2 = np.array(zc)[-1] + (np.array(zc)[-1]-np.array(zc)[-2])
				edges = np.concatenate((np.array([lim1]), edges))
				edges = np.concatenate((edges, np.array([lim2])))
				self.data.update({'bin_edges': edges})
				print 'Got bin centres.'
			if ids: 
				id = [i[0] for i in self.hdus[2].data[:nsamp]]
				self.data.update({'coadd_objects_id':np.array(id)})
				print 'Got object ids.'
			if modes: 
				md= [m[0] for m in self.hdus[3].data[:nsamp]]
				self.data.update({'pz_cond_mode':np.array(md)})
				print 'Got %d PDF modes.'%nsamp
			if means: 
				mn = [m[0] for m in self.hdus[4].data[:nsamp]]
				self.data.update({'pz_cond_mean':np.array(mn)})
				print 'Got %d PDF means.'%nsamp
			if pz_samples: 
				pz = [p[0] for p in self.hdus[5].data[:nsamp]]
				self.data.update({'z_samples':np.array(pz)})
				print 'Got %d redshift samples.'%nsamp
			if stdv: 
				std = [s[0] for s in self.hdus[6].data[:nsamp]]
				self.data.update({'pz_cond_std':np.array(std)})
				print 'Got PDF standard deviation for %d samples.'%nsamp

	def find_photoz_bin_edges(self, nbin, sel=None, binning='mode'):
		if sel==None:	sel = np.in1d(self.data['pz_cond_%s'%binning], self.data['pz_cond_%s'%binning])

		self.data['tomographic_bin_edges'] = [ min(self.data['pz_cond_%s'%binning][sel]) ]
		for i in range(1,nbin+1):
			edge = sp.percentile(self.data['pz_cond_%s'%binning][sel],i*100.0/nbin)
			self.data['tomographic_bin_edges'] += [edge]
		self.data['tomographic_bin_edges'] += [max(self.data['pz_cond_%s'%binning][sel])]
		self.data['tomographic_bin_edges'] = np.array(self.data['tomographic_bin_edges'])

	def evaluate_pz_tr(self, zmin, zmax, bin, selection_mask=None, binning='mode'):
		n_z = []

		print 'Using %s binning.'%binning
		if binning=='mode':
			binning = 'pz_cond_mean'
		elif binning=='mean':
			binning = 'pz_cond_mode'
		else:
			print 'Unrecognised binning option. Please choose binning = mean, mode'
			print 'Using default (mode).'

		# First apply the selection mask for the required subcatalogue
		if selection_mask==None:	selection_mask = np.in1d(self.data['pz_cond_mean'], self.data['pz_cond_mean'])

		# Then select samples from p(z_phot) within bin edges
		photoz_bin_filter = (self.data[binning]>zmin) & (self.data[binning]<zmax)
		for i in range(len(self.data['bin_edges'])-1):
			sel = selection_mask & photoz_bin_filter
			z1, z2 = self.data['bin_edges'][i], self.data['bin_edges'][i+1]
			# Third mask to define fine binning used to evaluate p(z_tr)
			sel = sel & (self.data['z_samples']>z1) & (self.data['z_samples']<z2)
			n_z += [len(self.data['z_samples'][sel])]
		n_z = np.array(n_z)

		# Normalise
		zf = np.linspace(self.data['bin_edges'][0], self.data['bin_edges'][-1], 2000)
		pzf = np.interp(zf, self.data['bin_centres'], n_z)
		p_z = pzf / sum(pzf * (zf[1]-zf[0]) )
		self.data.update({"z": zf})
		self.data.update({"pz_bin_%d"%bin: p_z})
		self.data.update({"nz_bin_%d"%bin: n_z})

	def output_pz(self, plot=True, save=False, save_dir=None):
		if save:		
			out={}
			import cPickle as pi
		for i in range(len(pz.data['tomographic_bin_edges'])+1):
			plt.plot(pz.data['z'], pz.data['pz_bin_%d'%i])
			if i!=0: plt.axvline(pz.data['tomographic_bin_edges'][i-1],color='k',linestyle='--')
			if i<len(pz.data['tomographic_bin_edges']): out['bin_%d'%i] = pz.data['pz_bin_%d'%i]
		out['z'] = pz.data['z']
		out['z_bin_edges'] = pz.data['tomographic_bin_edges']

	def stack_pdfs(self, zmin, zmax, nz, name):
		z = np.linspace(0,3,nz)

		sel = (self.data['point_z']>zmin) & (self.data['point_z']<zmax)
		point_z,sigma_z = self.data['point_z'][sel], self.data['sigma_z'][sel]

		p_z = np.exp(-1.0*(z-point_z[0])*(z-point_z[0])/2.0/sigma_z[0]/sigma_z[0])

		for i in range(1, len(point_z)):
			z_p, dz = point_z[i], sigma_z[i]
			# Define a Guassian PDF for each point redshift estimate 			
			p_z += np.exp(-1.0*(z-z_p)*(z-z_p)/2.0/dz/dz)

		# Normalise to 1 over all redshifts
		p_z = p_z / sum(p_z*(z[1]-z[0]))

		self.data.update({name : {'nz':p_z, 'z':z }})


class redmagic_calc(calculator):

	cosmology = {"h0":0.6727, "sigma8":0.82, "omega_m":0.3156, "omega_k":0.0, "w0":-1.0, "wa":0.0, "ns":0.9645, "omega_b":0.0491685}
	code="redmagic"
	def __init__(self, analytic=False, from_txt=False, filename=None, cat2=False, fhd="y1a1_gold_1.0.2b-full_redmapper_v6.4.11_redmagic_highdens_0.5-10.fit", fhl="y1a1_gold_1.0.2b-full_redmapper_v6.4.11_redmagic_highlum_1.0-04.fit"):
		import pyfits
		if filename and not from_txt:
			self.data = pyfits.getdata(filename)
			print "Loaded catalogue of %fM galaxies" %(len(self.data)/1e6)
			self.zspec = self.data["ZSPEC"]
			self.zphot = self.data["ZREDMAGIC"]
			self.zphot_sigma = self.data["ZREDMAGIC_E"]
			self.g, self.r, self.i, self.z = self.data["MABS"].T
		elif filename and from_txt:
			self.files=[filename]

		if cat2:
			self.highL=pyfits.getdata(fhl)
			self.highrho=pyfits.getdata(fhd)

		if analytic:
			print "Set up to generate mock redmagic n(z)"
			print "for details see Rozo et al 2015 arXiv:1507.05460v1.pdf"
			self.ncm = 10e-3
			self.zlim=(0.2,0.8)
			self.zphot = np.linspace(0.2, 0.8, 70000)
			self.zphot_sigma = 0.017 * (1.0+self.zphot)


	def redshift_to_cmdistance(self, z):
		"""Generic function to calculate the corresponding comoving distance bound
		for a given redshift limit."""
		c = constants["c"]
		zp = np.linspace(0.0, z, 200)
		integrand = cosmology["omega_m"]*(1.0+zp) * (1.0+zp) * (1.0+zp)
		integrand += cosmology["omega_k"]* (1.0+zp) * (1.0+zp)
		w0,wa = cosmology["w0"], cosmology["wa"]
		f_de = np.exp( -3.0*wa* (1.0-1/(1.0+zp)) ) * (1.0+zp)**(3.0*(1.0+w0+wa))
		integrand += (1.0-cosmology["omega_m"]-cosmology["omega_k"]) * f_de

		integrand = 1.0/np.sqrt(integrand)
		x = np.trapz(integrand, zp)

		return x


	def get_stacked_pdfs(self, ngal=-1, spec=False, edges=[0.,1.3]):
		if spec:
			sel = self.data["ZSPEC"]>0.
		else:
			 sel = np.ones_like(self.zphot).astype(bool)

		selbin = (self.zphot>edges[0]) & (self.zphot<edges[1])
		z_in_bin = self.zphot[selbin]
		# Define some phot-z samples within the limits of this bin
		if ngal!=-1:
			sel=np.random.randint(0,len(z_in_bin), ngal)
		x = np.linspace(0.0, 1.5, 200)
		n_z = np.zeros_like(x)		
		# Draw random point estimates (sampling from p(zphot))
		# Estimate a Gaussian p(ztr|zphot) for each one
		for i, z in enumerate(z_in_bin[sel]):
			n_z += gaus(x, z, self.zphot_sigma[selbin][sel][i])

		return x, n_z

	def get_equal_number_bins(self, nbin=3):
		edges=np.array([])
		for i in xrange(nbin+1):
			e = sp.percentile(self.zphot, i*100./nbin)
			edges = np.append(edges,e)
		return edges

	def get_n_of_z(self,nbin=3, integration_points=20000, bin_edges=None):
		if not bin_edges:
			print "Using equal number bins."
			bins = self.get_equal_number_bins(nbin=nbin)
		else:
			bins = bin_edges
		print "Calculated bin edges:", bins

		nz = []		
		for i in xrange(len(bins)-1):
			z, n = self.get_stacked_pdfs(ngal=integration_points, edges=[bins[i], bins[i+1]])
			nz += [n]
			print "Calculated bin %d"%i

		self.z = z
		self.edges=bins
		self.pz = nz

		return z, bins, nz

def save(filename, calculator):
	z = getattr(calculator, "z")
	pz = getattr(calculator, "pz")
	out=z
	for p in pz: 
		out=np.vstack((out,p))
	out = out.T
	np.savetxt(filename, out)

def plot(calculator, xmin=0.0, xmax=1.3, normalise=True, colour="purple", label=None, linestyle="-"):
		import pylab as plt
		pz = getattr(calculator, "pz")
		z = getattr(calculator, "z")
		for i, p in enumerate(pz):
			if normalise:
				norm = np.trapz(p, z)
			else:
				norm=1.0

			if label and (i==0):
				lab = label
			else:
				lab=None
			plt.plot(z, p/norm, colour, linestyle=linestyle, label=lab)

		plt.xlabel("redshift")
		plt.ylabel("$n(z)$")
		plt.xlim(xmin,xmax)




#	def calculate_ngal(self, zedges=None, area=1523.23244035):
#		NL=[]
#		area = area*60*60
#        for zmax in [2.6,2.5,2.3,2.18, 2.09,2.0,1.8,1.6,1.0]:
#        	f=zmax/2.0
#        	if not zedges:
#        		edges=np.linspace(0.3,.75*f,4)
#        		nbin=3
#        	else:
#        		edges=[zedges[0]] + list(np.array(zedges[1:])*f)
#        		nbin = len(edges)-1
#
#        	print "bin edges:",edges
#        	N=[]
#        	for i in range(nbin-1):
#        		sel = (self.highrho["zredmagic"]>edges[i]) &  (self.highrho["zredmagic"]<edges[i+1])
#        		N+=[len(self.highrho[sel])]
#        	for i in range(nbin-1,nbin):
#        		sel = (self.highL["zredmagic"]>edges[i]) &  (self.highL["zredmagic"]<edges[i+1])
#        		N+=[len(self.highL[sel])]
#        		NL+=[N]
                #print NL
        #return np.array(NL)/area



def gaus(x,x0,sig):
	return np.exp(-1.*(x-x0)*(x-x0)/2/sig/sig)/np.sqrt(2.*np.pi)/sig

def match_DES_sv_photoz(catfile, svfile="/home/samuroff/des_sv_wl_info.fits"):
	nz = pyfits.getdata(svfile)
	cat = pyfits.getdata(catfile)
	col=np.vstack((nz['ra'], nz['dec']))
	svcol=np.vstack((cat['ALPHAWIN_J2000'], cat['DELTAWIN_J2000']))
	tree=sps.KDTree(col.T)
	nn=tree.query(svcol.T)
	svind=nn[1]
	selected_sv=nz[svind]

	# Match to ~1 pixel in both ra and dec
	diffsel=(abs(selected_sv['ra']-cat['ALPHAWIN_J2000'])<7.5e-5) & (abs(selected_sv['dec']-cat['DELTAWIN_J2000'])<7.5e-5)
	zdes=selected_sv[diffsel]['mean_photoz']

	return selected_sv

class sv_pz(calculator):
	bin_info={
	    "skynet" : (np.linspace(0.005, 1.8, 201)),
	    "skynet_biased" : (np.linspace(0.005, 1.8, 201)),
	    "annz" : (np.linspace(0.00, 1.8, 181)),
	    "bpz" : (np.linspace(0.005, 2.505, 251)),
	    "raw_bpz" : (np.linspace(0.055, 2.555, 251)),
	    "tpz" : (np.linspace(0.005-(0.007474/2), 2.0-(0.07474/2), 201)),
	    "uberz" : (np.arange(0,4.0+0.01/2,0.01)),
	    }

	def __init__(self, filename, code):
		import glob
		self.files = glob.glob(filename)
		self.code = code.lower()

		self.edges=self.bin_info[self.code]
		self.z = (self.edges[1:] + self.edges[:-1])/2.0

		print "%d file(s) to load"%len(self.files)

	def extract_info(self, store=True):
		import fitsio
		dat =fitsio.FITS(self.files[0])

		self.columns = dat[1].get_colnames()
		self.npdf=0
		for c in self.columns:
			if "pdf" in c.lower():
				self.npdf+=1

		if store:
			self.data = dat


	def extract_coadd_ids(self):
		self.coadd_ids = self.dat[1]["coadd_objects_id"]
		self.coadd_ids = self.coadd_ids.read()

	def build_n_of_z(self, nbin,  bins=None, equal_number_bins=True, zlim=(0.0,2.5), frac=1.0, store=True):

		if hasattr(self,"zmean"):
			zmean = self.zmean
		else:
			zmean = self.data[1]["z_mean"].read()
			self.mask = np.isfinite(zmean)
			zmean = zmean[self.mask]
			self.zmean = zmean

		# Determine the tomographic bin edges
		if equal_number_bins:
			bins = self.get_equal_number_bins(nbin, zlim)
		else:
			bins = bins
			if bins == None:
				print "Error: if equal_number_bins!=T then bin edges must be specified."

		self.tomographic_bin_edges=bins

		print "Constructing %d tomographic bins with edges:"%nbin
		print bins

		self.pz = []

		# Then sum the PDFs within the bin limits to compute p(z) in each fine z bin
		# Loop over tomographic bins
		for i in xrange(len(bins)-1):
			b0,b1 = bins[i], bins[i+1]

			pz = np.zeros_like(self.z)
			print "Stacking PDFs in %d fine bins"%self.npdf

			# Then over the fine redshift bins
			for j in xrange(self.npdf):
				pdf_name = "pdf_%d"%j
				print j

				# Load the PDF values (or some subset of them) at the fine bin centre and add them together 
				nmax = len(zmean)
				if frac!=1.0:
					nmax *= frac
					nmax = int(nmax)

				if not hasattr(self, "pdf_%d"%j):
					p = self.data[1][pdf_name].read()[self.mask][:nmax]
					setattr(self, "pdf_%d"%j, p)
				else:
					p = getattr(self, "pdf_%d"%j)

				sel = (zmean[:nmax]>b0) & (zmean[:nmax]<b1)

				pz[j] = sum(p[sel])

			self.pz.append(pz)

	def calibrate_bpz(self):
		dz = 0.050
		if self.code=="bpz":
			print "Applying BPZ calibration dz=%3.2f"%dz
			self.z+=dz


	def get_equal_number_bins(self, nbin, zlim):
		if hasattr(self, "zmean"):
			zmean = self.zmean
		else:
			zmean = self.data[1]["z_mean"].read()
			self.zmean=zmean
		sel = (zmean>zlim[0]) & (zmean<zlim[1])

		return find_bin_edges(zmean[sel], nbin)

	def plot(self):
		import pylab as plt
		for p in self.pz:
			plt.plot(self.z, p, "purple")

		plt.xlabel("redshift")
		plt.ylabel("$n(z)$")

	def bias(self, deltaz):
		for i,p in enumerate(self.pz):
			print "Applying bias of %3.3f to bin %d"%(deltaz,i+1)
			self.pz[i] = self.translate(i,deltaz)

	def translate(self, bin, deltaz):
		"""Shift a redshift distribution upwards or downwards."""
		pz = self.pz[bin]
		zf = np.linspace(0, self.z.max(), 1000)

		pzfine=np.interp(zf, self.z, pz)
		resolution = (zf[1:]-zf[:-1])[0]
		shift_required = deltaz/resolution

		# Pad with zeros so there are no weird effects from nonzero end points
		# shifted from one end of the array to the other
		pzfine=[0.0]*100+list(pzfine)+[0.0]*100
		pzfineb=np.roll(np.array(pzfine), int(math.floor(0.5+shift_required)))
		pzfineb=pzfineb[100:-100]

		return np.interp(self.z,zf,pzfineb)

def find_bin_edges(x,nbins,w=None):
    """
    For an array x, returns the boundaries of nbins equal (possibly weighted by w) bins.
    Stolen from Troxel's code.
    """

    if w is None:
      xs=np.sort(x)
      r=np.linspace(0.,1.,nbins+1.)*(len(x)-1)
      return xs[r.astype(int)]

    fail=False
    ww=np.sum(w)/nbins
    i=np.argsort(x)
    k=np.linspace(0.,1.,nbins+1.)*(len(x)-1)
    k=k.astype(int)
    r=np.zeros((nbins+1))
    ist=0
    for j in xrange(1,nbins):
      # print k[j],r[j-1]
      if k[j]<r[j-1]:
        print 'Random weight approx. failed - attempting brute force approach'
        fail=True
        break
      w0=np.sum(w[i[ist:k[j]]])
      if w0<=ww:
        for l in xrange(k[j],len(x)):
          w0+=w[i[l]]
          if w0>ww:
            r[j]=x[i[l]]
            ist=l
            break
      else:
        for l in xrange(k[j],0,-1):
          w0-=w[i[l]]
          if w0<ww:
            r[j]=x[i[l]]
            ist=l
            break

    if fail:

      ist=np.zeros((nbins+1))
      ist[0]=0
      for j in xrange(1,nbins):
        wsum=0.
        for k in xrange(ist[j-1].astype(int),len(x)):
          wsum+=w[i[k]]
          if wsum>ww:
            r[j]=x[i[k-1]]
            ist[j]=k
            break

    r[0]=x[i[0]]
    r[-1]=x[i[-1]]

    return r


class combined_pz(calculator):
	def __init__(self, pzlist, bins):
		print "initialised calculator object to create hybrid n(z) using %d bins (%d unique input n(z))"%(len(pzlist), len(np.unique(pzlist)))
		self.pzs=pzlist
		self.codes=[]
		for p in pzlist:
			self.codes.append(p.code)

		print self.codes
		self.bins=bins

	def make(self):
		"""Put together the input p(z) in the bin/code combination specified."""
		self.pz=[]
		self.z=np.array([1.0e6])

		for p in self.pzs:
			if (p.z.max() < self.z.max()):
				self.z = p.z
				c = p.code

		print "Adopting redshift sampling of %s. Will interpolate the other p(z) to the same points."%c


		for i, input_pz in enumerate(self.pzs):
			interpolated_pz = np.interp(self.z, input_pz.z, input_pz.pz[self.bins[i]])
			self.pz.append(interpolated_pz)


class nofz:
	def __init__(self, filename):
		self.z = np.loadtxt(filename).T[0]

		self.bins=[]
		for nz in np.loadtxt(filename).T[1:]: self.bins.append(nz/np.trapz(nz,self.z))

	def generate_interpolators(self):
		print "Setting up p(z) interpolators"

		self.pz=[]

		for nz in self.bins:
			self.pz.append(spi.interp1d(self.z, nz, kind="cubic"))

	def assign_galaxies_to_bins(self, z_true, verbose=True, zlim=True):

		self.true_z = z_true
		if zlim:
			self.bin_allocation = np.zeros_like(z_true)
			z_true[z_true<self.z.min()] = self.z.min()
			upper_redshift_limit = (z_true<self.z.max())
			z_true = z_true[upper_redshift_limit] 


		# Evaluate each of the p(z) at a specific COSMOS (true) redshift
		if verbose:
			print "Evaluating p(z)"
		peval = [pev(z_true) for pev in self.pz]

		ptot = np.array(peval).sum(axis=0)

		# Divide through by the total weight at this true redshift
		# This should give us a probablility of finding each COSMOS galaxy in each 
		# of the redshift bins
		if verbose:
			print "Converting to bin probabilities"
		for i in xrange(len(peval)): peval[i] = peval[i]/ptot

		# Now generate a random number for each COSMOS galaxy
		random = np.random.rand(peval[0].size)

		# Finally map that onto the bins to assign a DES redshift bin
		if verbose:
			print "Allocating galaxies across %d bins"%len(peval)
		edges = [np.array(peval)[:i].sum(axis=0) for i in xrange(len(peval)+1)]

		ntot = peval[0].size
		bin_allocation=np.zeros_like(peval[0])-1
		for i, (lower,upper) in enumerate(zip(edges[:-1], edges[1:])):
			selection = (random<upper) & (random>lower)
			bin_allocation[selection] = (i+1)
			ngal = bin_allocation[selection].size
			if verbose:
				print "%d/%d (%3.3f percent ) galaxies assigned to bin %d"%(ngal,ntot,100.*ngal/ntot,i+1)

		self.bin_allocation[upper_redshift_limit] = bin_allocation

		return self.bin_allocation

	def histograms(self, filename="/home/samuroff/shear_pipeline/plot_dump/hoopoe/hoopoe_redshift_bin_allocation.png"):
		import pylab as plt
		labels=["True COSMOS Redshift", None,None,None]
		labels2=["DES BPZ", None,None,None]

		for i, n in enumerate(self.bins):
			plt.plot(self.z, n, lw=2.5, color=colors[i], alpha=0.5, label=labels2[i])

		for i in [0,1,2,3]:
			plt.hist(self.true_z[self.bin_allocation==i+1], lw=2.5, color=colors[i], normed=1, bins=60, histtype="step" , ls="dashed", label=labels[i] )

		plt.xlabel("Redshift $z$")
		plt.xlim(0,1.8)
		plt.legend()
		plt.yticks(visible=False)

		plt.savefig(filename)









		
