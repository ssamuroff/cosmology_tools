import numpy as np
import astropy.io.fits as pf
import fitsio as fio
import glob
#import treecorr as tc
import meds, fitsio
import galsim
import os, math
import py3shape.i3meds as i3meds
import py3shape.utils as utils
import py3shape.options as i3opt
import tools.diagnostics as di
from scipy.interpolate import Rbf
import tools.arrays as arr
import tools.neighbours as nf
from plots import im3shape_results_plots as i3s_plots
import py3shape as p3s


BLACKLIST_CACHE = {}
PSFEX_CACHE = {}
PSFEX_CACHE_FILENAME = ""

names = {"snr": "snr", "rgp" : "mean_rgpp_rp", "e1": "e1", "e2" : "e2", "iterations":"levmar_iterations", "psf_e1" : "mean_psf_e1", "psf_e2" : "mean_psf_e2", "psf_fwhm": "mean_psf_fwhm", "nflag": "neighbour_flag", "dphi":"dphi"}

sersic_indices={"disc": 1, "bulge": 4}

class shapecat(i3s_plots):
	"""Class for handling im3shape results and,
	   if relevant, truth catalogues.
	"""
	def __init__(self, res=None, truth=None, coadd=False, fit="disc", noisefree=False, res_arr=None):
		if res is None and truth is None:
			print "Please specify at least one of res=, truth="
		self.res_path = res
		self.truth_path = truth
		self.coadd_path= coadd

		self.fit = fit
		self.noisefree=noisefree

		self.corr={}

		if res_arr is not None:
			self.res=res_arr

		print "Initialised."

	def load_from_array(self, array, name="res"):
		setattr(self, name, array)

	def load_from_cat(self, cat, name="res"):
		if name is not "all":
			table = getattr(cat, name) 
			setattr(self, name, table)
			setattr(self, "%s_path"%name, getattr(cat, "%s_path"%name))
			self.files = cat.files
			try:
				self.indices = cat.indices
			except: 
					print "No indices column"

		else:
			if hasattr(cat, "truth") and name=="all":
				self.truth = cat.truth
				try:
					self.truth1 = cat.truth1
					self.truth2 = cat.truth2
				except:
					print "Once truth catalogue found (no splits)"

			if hasattr(cat, "res") and name=="all":
				self.res = cat.res
				try:
					self.res1 = cat.res1
					self.res2 = cat.res2
				except:
					print "Once results catalogue found (no splits)"


	def do_rbf_interpolation(self, bias, cat):

		# If the arrays here  are sufficiently large we may need to apply the interpolation to
		# the catalogue in parts
		if cat.size>25e6:
			nslice = cat.size/2
			if bias=="m":
				return np.hstack(( self.rbf_interp_m(np.log10(cat["snr"][:nslice])/self.fx, np.log10(cat["mean_rgpp_rp"][:nslice])/self.fy), self.rbf_interp_m(np.log10(cat["snr"][nslice:])/self.fx, np.log10(cat["mean_rgpp_rp"][nslice:])/self.fy) ))
			elif bias=="a":
				return np.hstack(( self.rbf_interp_a(np.log10(cat["snr"][:nslice])/self.fx, np.log10(cat["mean_rgpp_rp"][:nslice])/self.fy), self.rbf_interp_a(np.log10(cat["snr"][nslice:])/self.fx, np.log10(cat["mean_rgpp_rp"][nslice:])/self.fy) ))
		else:
			if bias=="m":
				return self.rbf_interp_m(np.log10(cat["snr"])/self.fx, np.log10(cat["mean_rgpp_rp"])/self.fy)
			elif bias=="a":
				return self.rbf_interp_a(np.log10(cat["snr"])/self.fx, np.log10(cat["mean_rgpp_rp"])/self.fy)

	def get_weights_grid(self, sbins=10, rbins=5,rlim=(1,3), slim=(10,1000), binning="equal_number"):
		print 'calculating weights on a grid'

		data = self.res
		

		sel_lim = (data["snr"]>slim[0]) & (data["snr"]<slim[1]) & (data["mean_rgpp_rp"]>rlim[0]) & (data["mean_rgpp_rp"]<rlim[1])
		data = data[sel_lim]

		if isinstance(binning,str) : 
			if binning.lower()=="uniform":
				snr_edges = np.logspace(np.log10(slim[0]),np.log10(slim[1]),sbins+1)
				rgp_edges = np.linspace(rlim[0],rlim[1],rbins+1)
			elif binning.lower()=="equal_number":
				snr_edges = di.find_bin_edges(np.log10(data["snr"]), sbins)
				rgp_edges = di.find_bin_edges(np.log10(data["mean_rgpp_rp"]), rbins)

		
		snr_centres = (snr_edges[1:]+snr_edges[:-1])/2.0
		rgp_centres = (rgp_edges[1:]+rgp_edges[:-1])/2.0

		w_list=[]
		w_grid=[]

		print "Will do dynamic binning in SNR"

		for i in xrange(len(rgp_edges)-1):
			snr_samp = data["snr"][(np.log10(data['mean_rgpp_rp']) > rgp_edges[i]) & (np.log10(data['mean_rgpp_rp']) < rgp_edges[i+1])]
			snr_edges=di.find_bin_edges(np.log10(snr_samp), sbins)
			for j in xrange(len(snr_edges)-1):
				empty=False
				print "bin %d %d  snr = [%2.3f-%2.3f] rgpp/rp = [%2.3f-%2.3f]"%(j, i, 10**snr_edges[j], 10**snr_edges[j+1], 10**rgp_edges[i], 10**rgp_edges[i+1] )

				# Select in bins of snr and size
				select = (np.log10(data['snr']) > snr_edges[j]) & (np.log10(data['snr']) < snr_edges[j+1]) & (np.log10(data['mean_rgpp_rp']) > rgp_edges[i]) & (np.log10(data['mean_rgpp_rp']) < rgp_edges[i+1])
				ngal = np.nonzero(select.astype(int))[0].size

				# Raise an error if there are too few simulated galaxies in a given bin
				if ngal < 60:
					print "Warning: <100 galaxies in bin %d, %d (ngal=%d)"%(i,j, ngal)
					empty=False
				if ngal==0:
					print "Warning: no galaxies in bin %d, %d "%(i,j)
					empty=True

				vrgp_mid = rgp_centres[i]
				vsnr_mid = snr_centres[j]
				vrgp_min = rgp_edges[i]
				vsnr_min = snr_edges[j]
				vrgp_max = rgp_edges[i+1]
				vsnr_max = snr_edges[j+1]

				if ngal==0:
					print "Warning: no galaxies in bin %d, %d"%(i,j)
					w_list.append([i,j,ngal, 10**vrgp_min, 10**vrgp_max, 10**vsnr_min, 10**vsnr_max, 0])
					continue

				w1 = di.compute_weight(data["e1"][select], verbose=False)
				w2 = di.compute_weight(data["e2"][select], verbose=False)
				w = (w1+w2)/2.0
				
				filename_str = 'snr%2.2f.rgpp%2.2f' % (10**vsnr_mid,10**vrgp_mid)
				w_list.append([j, i, ngal, 10**vrgp_min, 10**vrgp_max, 10**vsnr_min, 10**vsnr_max, w])

		dt = [("j", int), ("i", int), ("ngal", int), ("rgp_lower", float), ("rgp_upper", float), ("snr_lower", float), ("snr_upper", float), ("weight", float),]
		arr_weights = np.zeros(np.array(w_list).T[0].size, dtype=dt)

		for a,b in enumerate(arr_weights.dtype.names):
			arr_weights[b] = np.array(w_list).T[a]

		filename_table_bias = '/home/samuroff/weights_table.fits'
		
		import fitsio as fi

		out = fi.FITS(filename_table_bias, "rw")
		out.write(arr_weights)
		out[-1].write_key("EXTNAME", "im3shapev1_weights")
		out.close()
		print 'saved %s'%filename_table_bias

		self.weights_table = arr_weights

	def interpolate_weights_grid(self, nslice=100):
		sel = (np.isfinite(np.log10(self.res["snr"])) & np.isfinite(np.log10(self.res["rgpp_rp"])))
		x = np.sqrt(self.weights_table["snr_lower"]*self.weights_table["snr_upper"])
		y = np.sqrt(self.weights_table["rgp_lower"]*self.weights_table["rgp_upper"])
		fx = np.log10(self.res["snr"][sel]).max()
		fy = np.log10(self.res["rgpp_rp"][sel]).max()
		ngal = self.res["rgpp_rp"][sel].size
		interpolator = Rbf(np.log10(x)/fx, np.log10(y)/fy, self.weights_table["weight"], smooth=3, function="multiquadric")

		print "Will do interpolation in %d sections."%nslice

		galaxy_weights=np.zeros(ngal)
		start=0
		end=ngal/nslice

		for i in xrange(nslice):
			print i, start, end
			galaxy_weights[start:end] = interpolator(np.log10(self.res["snr"][sel][start:end])/fx, np.log10(self.res["rgpp_rp"][sel][start:end])/fy)
			start += ngal/nslice
			end += ngal/nslice 
			
		print "Calculated weights for %s galaxies"%galaxy_weights.size
		print "Mean weight : %2.4f, Min : %2.4f, Max : %2.4f"%(galaxy_weights.mean(), galaxy_weights.min(), galaxy_weights.max())

		return sel, galaxy_weights


	def load(self, res=True, truth=False, epoch=False, coadd=False, prune=False, cols=[None,None], postprocessed=True, keyword="DES", apply_infocuts=True, ext=".fits", match=[], ntiles=None):
		
		if res and (not hasattr(self, "res")):
			if "%s"%ext in self.res_path:
				files=[self.res_path]
			else:
				files = glob.glob("%s/*%s"%(self.res_path,ext))

			if ntiles is not None:
				files = files[:ntiles]
			if len(match)>0:
				tiles = [os.path.basename(f)[:12] for f in files]
				select=np.array([(t in match )for t in tiles])
				files=np.array(files)[select]
				ntiles = files.size

			print "%s/*%s"%(self.res_path,ext)
			single_file=False
			print "loading %d results file(s) from %s"%(len(files),self.res_path)

			if self.noisefree and apply_infocuts:
				self.res = pf.getdata(files[0])
				tmp, noise_free_infocuts = self.get_infocuts(exclude=["chi"], return_string=True)

				#noise_free_infocuts = noise_free_infocuts.replace("cuts= ((", "cuts= ((%s['chi2_pixel']>0.004) & (%s['chi2_pixel']<0.2) & (")

			if len(files)>1:
				if apply_infocuts and self.noisefree:
					self.res, self.files, i = di.load_results(res_path =self.res_path, format=ext[1:], cols=cols[0], apply_infocuts=False, additional_cuts=noise_free_infocuts, keyword=keyword, postprocessed=postprocessed, return_filelist=True, match=match)
				else:
					self.res, self.files, i = di.load_results(res_path =self.res_path, format=ext[1:], cols=cols[0], ntot=len(files) ,apply_infocuts=apply_infocuts, keyword=keyword, postprocessed=postprocessed, return_filelist=True, match=match)
			else:
				if ext.lower()==".fits":
					self.res = fio.FITS(files[0])[1].read(columns=cols[0])
				elif ext.lower()==".txt":
					self.res = np.genfromtxt(files[0], names=True)
				self.files=files
				single_file=True
				i=None

			self.indices=i


		if truth:
			if ".fits" in self.truth_path:
				files=[self.truth_path]
			else:
				files = glob.glob("%s/*.fits*"%self.truth_path)
			single_file=False

			print "loading truth files from %s"%self.truth_path
			if len(files)>1:
				if res:
					self.truth = di.load_truth(truth_path=self.truth_path, cols=cols[1], apply_infocuts=apply_infocuts, match=self.files, ind=self.indices, res=self.res)
				else:
					self.truth = di.load_truth(truth_path=self.truth_path, cols=cols[1], apply_infocuts=apply_infocuts)
			else:
				self.truth = pf.getdata(files[0])

		if truth and res:
			self.res, self.truth = di.match_results(self.res,self.truth, name1="DES_id", name2="coadd_objects_id", unique=False)
			if ("ra" in self.res.dtype.names): 
				if not (self.res["ra"]==self.truth["ra"]).all():
					self.res["ra"] = self.truth["ra"]
					self.res["dec"] = self.truth["dec"]

			print "Found catalogue of %d objects after matching to truth table"%len(self.res)

		if coadd:
			if ".fits" in self.coadd_path:
				files=[self.coadd_path]
			else:
				files = glob.glob("%s/*cat*.fits"%self.coadd_path)
			single_file=False

			print "loading coadd catalogue files from %s"%self.coadd_path
			if len(files)>1:
				print "update code..."
			else:
				self.coadd = pf.getdata(files[0])

			ids = self.res["row_id"]-1
			self.coadd= self.coadd[ids]

		if epoch:
			path = self.res_path.replace("main", "epoch")
			try:
				self.epoch = di.load_epoch(path)
			except:
				self.epoch = di.load_epoch(path.replace("bord", "disc"))



#		if hasattr(self, "truth"):
#			sel = self.truth["sextractor_pixel_offset"]<1.0
#			self.truth = self.truth[sel]
#			if hasattr(self, "res"):
#				self.res = self.res[sel]
#			if coadd:
	#			self.coadd = self.coadd[sel]

	def add_bpz_cols(self, fil="/share/des/disc3/samuroff/y1/photoz/bpz/NSEVILLA_PHOTOZ_TPL_Y1G103_1_bpz_highzopt_2_9_16.fits", array=None, exclude=None):
		print "Loading BPZ results from %s"%fil
		if array is None:
			bpz = fio.FITS(fil)[1].read()
		else:
			bpz = array

		self.res, bpz = di.match_results(self.res, bpz, name1="coadd_objects_id", name2="coadd_objects_id")

		for colname in bpz.dtype.names:
			if colname=="coadd_objects_id":
				continue
			if exclude is not None:
				if colname in exclude:
					continue
			else:
				print "Adding column: %s"%colname 
				self.res = add_col(self.res, colname, bpz[colname])

	def get_bulge_fraction(self, fcat, split_type="bpz", return_mask=False):
		if "bpz" in split_type:
			sel = fcat["T_B"]<2.0
		if "im3shape" in split_type:
			sel = fcat["bulge_flux"]!=0.0
		if return_mask:
			return 1.0*fcat[sel].shape[0]/fcat.shape[0], sel
		else:
			return 1.0*fcat[sel].shape[0]/fcat.shape[0]

	def get_positions(self, postype, selection_mask=None):
		"""Extract the position coordinates of objects, either from a truth
		   table or from an im3shape results table."""
		if postype not in ["world", "pixel"]:
			"Please choose a coordinate system ('world' or 'pixel')"

		if selection_mask is None:
			selection_mask = np.ones(self.res["coadd_objects_id"].size).astype(bool)

		if postype=="world":
			xname,yname = "ra","dec"
		elif postype=="pixel":
			xname,yname = "x_image","y_image"

		print "Obtained position data from",

		try:
			x,y = self.truth[xname][selection_mask],self.truth[yname][selection_mask]
			print "truth table"
		except:
			try:
				x,y = self.res[xname][selection_mask],self.res[yname][selection_mask]
				print "results table"
			except:
				x,y = self.coadd[xname.upper()][selection_mask],self.coadd[yname.upper()][selection_mask]
				print "coadd catalogue"

		return x,y


	def hist(self, param, normed=1, histtype="step", alpha=1.0, truth=False, res=True, nbins=25, xlim_upper=None, xlim_lower=None, label=None, colour="purple", infoflags=False, errorflags=False, neighbourflags=False, starflags=False, linestyle="solid", extra_cuts=None):
		"""Plot out a histogram of a particular parameter."""
		colname=names[param]
		sel= self.get_mask(info=infoflags, error=errorflags, stars=starflags, neighbours=neighbourflags)
		try:
			bins = np.linspace(xlim_lower, xlim_upper, nbins)
			plt.xlim(xlim_lower,xlim_upper)
		except:
			bins=nbins

		if not truth:
			data=self.res
		else:
			data=self.truth

		if extra_cuts is None:
			extra_cuts = np.ones_like(data)
		else:
			sel = sel & extra_cuts
		plt.hist(data[colname][sel], bins=bins, histtype=histtype, color=colour, normed=normed, lw=2.0, label=label, linestyle=linestyle, alpha=alpha)

	def get_mask(self, info=True, error=False, stars=False, neighbours=False):
		"""Get the selection mask for either error or info flag cuts."""
		if info:
			seli = self.res["info_flag"]==0
		else:
			seli = np.ones_like(self.res["e1"]).astype(bool)

		if error:
			sele = self.res["error_flag"]==0
		else:
			sele = np.ones_like(self.res["e1"]).astype(bool)

		if stars:
			sels = self.truth["star_flag"]==0
		else:
			sels = np.ones_like(self.res["e1"]).astype(bool)

		if neighbours:
			seln = self.truth["neighbour_flag"]==0
		else:
			seln = np.ones_like(self.res["e1"]).astype(bool)

		sel = (seli & sele & sels & seln)

		print "Mask leaves %d/%d objects."%(len(self.res[sel]), len(self.res))

		return sel

	def get_infocuts(self, exclude="None", return_string=False):
		cutstring = "cuts= ("
		count=0
		if exclude is None:
			exclude=[]
		for i, name in enumerate(info_cuts):
			exc=False
			for ex in exclude:
				if ex in name:
					exc=True
			if not exc:
				if count>0:
					cutstring+= "& %s "%name
				else:
					cutstring+= "%s "%name
				print i, name
				count+=1
		cutstring+=")"
		print cutstring
		exec cutstring

		if not return_string:
			return cuts

		else:
			return cuts, cutstring.replace("self.res", "%s")

	def rotate(self, e1_new, e2_new):
		"""Rotate the ellipticities into a reference frame defined by another set of shapes. 
		   Generally this will probably be the PSF or nearest neighbour frame.
		"""

		e = e1_new + 1j*e2_new
		# Actually this is 2phi
		phi = np.angle(e)

		self.res["e1"] = self.res["e1"]*np.cos(phi) + self.res["e2"]*np.sin(phi)
		self.res["e2"] = -1*self.res["e1"]*np.sin(phi) + self.res["e2"]*np.cos(phi)
		self.res["mean_psf_e1_sky"] = self.res["mean_psf_e1_sky"]*np.cos(phi) + self.res["mean_psf_e2_sky"]*np.sin(phi)
		self.res["mean_psf_e2_sky"] = -1*self.res["mean_psf_e1_sky"]*np.sin(phi) + self.res["mean_psf_e2_sky"]*np.cos(phi)
		self.truth["intrinsic_e1"] = self.truth["intrinsic_e1"]*np.cos(phi) + self.truth["intrinsic_e2"]*np.sin(phi)
		self.truth["intrinsic_e2"] = -1*self.truth["intrinsic_e1"]*np.sin(phi) + self.truth["intrinsic_e2"]*np.cos(phi)
		self.truth["true_g1"] = self.truth["true_g1"]*np.cos(phi) + self.truth["true_g2"]*np.sin(phi)
		self.truth["true_g2"] = -1*self.truth["true_g1"]*np.sin(phi) + self.truth["true_g2"]*np.cos(phi)

	def apply_infocuts(self):
		n0 = len(self.res)
		inf = self.res["info_flag"]==0
		self.res = self.res[inf]
		if hasattr(self,"truth"):
			self.truth=self.truth[inf]
		if hasattr(self,"coadd"):
			self.coadd=self.coadd[inf]
		print "%d/%d (%f percent) survived info cuts."%(len(self.res), n0, 100.*len(self.res)/n0)
		return 0

	def tune_neighbour_pool(self, task_number=0, tiles=[]):
		self.ncat.tune(self, task_number=task_number, tiles=tiles)

	def fast_neighbour_search(self, truth=None, repetition=False):
		if truth is None:
			truth = self.truth_path

		self.tiles = np.unique(self.res["tilename"]).astype("S12")
		object_tiles = self.res["tilename"].astype("S12")

		columns=['id', 'DES_id', 'cosmos_ident', 'cosmos_photoz', 'spectral_type', 'flags', 'star_flag', 'nearest_neighbour_pixel_dist', 'nearest_neighbour_index', 'true_g1', 'true_g2', 'intrinsic_e1', 'intrinsic_e2', 'intrinsic_e1_hsm', 'intrinsic_e2_hsm', 'hlr',  'ra', 'dec', 'mag', 'flux', 'nexp', 'mean_psf_e1', 'mean_psf_e2', 'mean_psf_fwhm']

		if not repetition:
			neighbours=[]
		else:
			neighbours=np.zeros_like(self.truth)

		print "Searching for neighbours in %d truth tables."%self.tiles.size

		for i, tile in enumerate(self.tiles):
			
			select = (object_tiles==tile)
			print i+1, tile, object_tiles[select].size

			filename = glob.glob("%s/%s*.fz"%(truth,tile))
			if len(filename)==0:
				print "Truth table is missing."
				continue
			if len(filename)>1:
				print "Warning - multiple truth tables found." 
				print "Will use %s"%filename[0]

			filename = filename[0]

			truth_table = fitsio.FITS(filename)[1].read(columns=columns)
			required = self.truth["nearest_neighbour_index"][select] 
	
			if repetition:
				finder=np.arange(0,len(truth_table), 1)
				indices=np.array([ finder[truth_table["DES_id"]==ni][0] for ni in required])
				neighbours[select] = truth_table[indices]
			else:
				extract = np.in1d(truth_table["DES_id"],required)  
				neighbours.append(truth_table[extract])

		if not repetition:
			neighbours = np.concatenate(neighbours)

		self.ncat = nf.neighbours(truth=neighbours)

		return neighbours


	def get_neighbours(self):
		import copy
		from sklearn.neighbors import NearestNeighbors
		import scipy.spatial as sps

		fulltruth = di.load_truth(self.truth_path)

		import fitsio as fi
		reference=fitsio.FITS("/share/des/disc6/samuroff/y1/hoopoe/y1a1-v2.2_10/meds/y1a1_positions.fits")[1].read()
		fulltruth,ref = di.match_results(fulltruth,reference, name1="coadd_objects_id", name2="DES_id")
		fulltruth["ra"]=ref["ra"]
		fulltruth["dec"]=ref["dec"]
		self.truth,ref = di.match_results(self.truth,reference, name1="coadd_objects_id", name2="DES_id")
		self.truth["ra"]=ref["ra"]
		self.truth["dec"]=ref["dec"]
		self.truth,self.res = di.match_results(self.truth,self.res, name1="coadd_objects_id", name2="DES_id")


		meds_path=self.truth_path.replace("truth", "meds/*/*")
		meds_info = di.get_pixel_cols(meds_path)
		pool_of_possible_neighbours,fulltruth = di.match_results(meds_info,fulltruth, name1="DES_id", name2="coadd_objects_id" )
		fulltruth = arr.add_col(fulltruth, "ix", pool_of_possible_neighbours["ix"])
		fulltruth = arr.add_col(fulltruth, "iy", pool_of_possible_neighbours["iy"])
		try:
			fulltruth = arr.add_col(fulltruth, "tile", pool_of_possible_neighbours["tile"])
		except:
			fulltruth["tile"] = pool_of_possible_neighbours["tile"]

		objects_needing_neighbours,self.truth = di.match_results(meds_info,self.truth, name1="DES_id", name2="coadd_objects_id" )
		self.truth = arr.add_col(self.truth, "ix", objects_needing_neighbours["ix"])
		self.truth = arr.add_col(self.truth, "iy", objects_needing_neighbours["iy"])
		try:
			self.truth = arr.add_col(self.truth, "tile", objects_needing_neighbours["tile"])
		except:
			self.truth["tile"] = objects_needing_neighbours["tile"]



		cut=(fulltruth["sextractor_pixel_offset"]<1.0) & (fulltruth["ra"]!=0.0)
		fulltruth = fulltruth[cut]
		pool_of_possible_neighbours = pool_of_possible_neighbours[cut]

		indices = np.zeros(self.res.size)
		distances = np.zeros(self.res.size)
		lookup = np.linspace(0,fulltruth.size-1, fulltruth.size).astype(int)

		tiles = np.unique(self.truth["tile"]) 

		for it in tiles:
			print "Matching in pixel coordinates, tile %s"%it
			sel0 = pool_of_possible_neighbours["tile"]==it
			sel1 = objects_needing_neighbours["tile"]==it

			# All positions where an object was simulated
			# Restrict the search to this tile
			x_pool = fulltruth["ra"][sel0] #pool_of_possible_neighbours["ix"][sel0]
			y_pool = fulltruth["dec"][sel0] #pool_of_possible_neighbours["iy"][sel0]
			xy_pool=np.vstack((x_pool,y_pool))

			# Positions of those objects for which we have im3shape results
			# We want to find neighbours for these objects
			x_tar = self.truth["ra"][sel1]
			y_tar = self.truth["dec"][sel1]
			xy_tar=np.vstack((x_tar,y_tar))

			# Build a tree using the pool
			nbrs = NearestNeighbors(n_neighbors=2, algorithm='kd_tree', metric="euclidean").fit(xy_pool.T)
			# Query it for the target catalogue
			d,i = nbrs.kneighbors(xy_tar.T)
			distances[sel1], indices[sel1] = d.T[1], lookup[sel0][i.T[1]]

		neighbour_cat = copy.deepcopy(self)

		neighbour_cat.res["id"]= fulltruth[indices.astype(int)]["DES_id"]
		neighbour_cat.res["coadd_objects_id"]= fulltruth[indices.astype(int)]["DES_id"]
		neighbour_cat.res["e1"]= fulltruth[indices.astype(int)]["intrinsic_e1"]+fulltruth[indices.astype(int)]["true_g1"] 
		neighbour_cat.res["e2"]= fulltruth[indices.astype(int)]["intrinsic_e2"]+fulltruth[indices.astype(int)]["true_g2"]
		np.putmask(neighbour_cat.res["e1"], neighbour_cat.res["e1"]<-1, fulltruth[indices.astype(int)]["mean_psf_e1"])
		np.putmask(neighbour_cat.res["e2"], neighbour_cat.res["e2"]<-1, fulltruth[indices.astype(int)]["mean_psf_e2"])
		neighbour_cat.res["ra"]= fulltruth[indices.astype(int)]["ra"]
		neighbour_cat.res["dec"]= fulltruth[indices.astype(int)]["dec"]
		neighbour_cat.truth= fulltruth[indices.astype(int)]
		neighbour_cat.truth["nearest_neighbour_pixel_dist"] = distances

		return neighbour_cat

	def match_to_faint(self):
		import copy
		from sklearn.neighbors import NearestNeighbors
		import scipy.spatial as sps

		faint = di.load_truth(self.truth_path, faint=True, add_tilename_col=True)

		indices = np.zeros(self.res.size)
		distances = np.zeros(self.res.size)
		lookup = np.linspace(0,faint.size-1, faint.size).astype(int)

		tiles = np.unique(self.truth["tilename"])

		for it in tiles:
			print "Matching in pixel coordinates, tile %s"%it
			sel0 = faint["tilename"]==it
			sel1 = self.truth["tilename"]==it
			

			# All positions where an object was simulated
			# Restrict the search to this tile
			x_pool = faint["ra"][sel0] #pool_of_possible_neighbours["ix"][sel0]
			y_pool = faint["dec"][sel0] #pool_of_possible_neighbours["iy"][sel0]
			xy_pool=np.vstack((x_pool,y_pool))

			# Positions of those objects for which we have im3shape results
			# We want to find neighbours for these objects
			x_tar = self.truth["ra"][sel1]
			y_tar = self.truth["dec"][sel1]
			xy_tar=np.vstack((x_tar,y_tar))

			# Build a tree using the pool
			nbrs = NearestNeighbors(n_neighbors=2, algorithm='kd_tree', metric="euclidean").fit(xy_pool.T)
			# Query it for the target catalogue
			d,i = nbrs.kneighbors(xy_tar.T)
			distances[sel1], indices[sel1] = d.T[1], lookup[sel0][i.T[1]]

		neighbour_cat = copy.deepcopy(self)

		neighbour_cat.res["id"]= faint[indices.astype(int)]["associated_object"]
		neighbour_cat.res["coadd_objects_id"]= faint[indices.astype(int)]["associated_object"]
		neighbour_cat.res["e1"]= faint[indices.astype(int)]["intrinsic_e1"]+faint[indices.astype(int)]["true_g1"] 
		neighbour_cat.res["e2"]= faint[indices.astype(int)]["intrinsic_e2"]+faint[indices.astype(int)]["true_g2"]
		neighbour_cat.res["ra"]= faint[indices.astype(int)]["ra"]
		neighbour_cat.res["dec"]= faint[indices.astype(int)]["dec"]
		neighbour_cat.truth= faint[indices.astype(int)]
		neighbour_cat.truth = arr.add_col(neighbour_cat.truth, "nearest_neighbour_pixel_dist", distances)

		return neighbour_cat

	def get_tangential_shear(self, x0, y0, x, y):
		dx = x-x0
		dy = y-y0
		sel = np.invert((dy==0) | (dx==0))

		theta = np.arctan2(dy[sel],dx[sel])

		et = -1*(self.res["e1"][sel]*np.cos(2*theta) + self.res["e2"][sel]*np.sin(2*theta))
		ex = -1*(self.res["e2"][sel]*np.cos(2*theta) - self.res["e1"][sel]*np.sin(2*theta))

		return et, ex, sel

	def get_phi_col(self,neighbour_cat):
		print "Computing ellipticity-ellipticity misalignment between each object and its nearest neighbour."
		eres = self.res["e1"]+ 1j*self.res["e2"]
		try:
			eneigh = neighbour_cat.res["e1"] + 1j*neighbour_cat.res["e2"]
		except:
			eneigh = neighbour_cat.truth["intrinsic_e1"] + 1j*neighbour_cat.truth["intrinsic_e2"]

		phi_res = np.angle(eres)
		phi_neigh = np.angle(eneigh)
		dphi = (phi_res - phi_neigh)/2

		# Enforce limits at +-pi, then +-pi/2
		#sel1=dphi>np.pi
		#sel2=dphi<-1.0*np.pi
#
		#dphi[sel1] = -1.0*np.pi + (dphi[sel1]-np.pi)
		#dphi[sel2] = 1.0*np.pi + (dphi[sel2]+np.pi)

		sel1=dphi>np.pi/2
		sel2=dphi<-1.0*np.pi/2
		
		dphi[sel1] = np.pi/2 - (dphi[sel1] - np.pi/2)
		dphi[sel2] = -1.*np.pi/2 - (dphi[sel2] + np.pi/2)

		dphi/=np.pi

		self.res = arr.add_col(self.res,"dphi",dphi)

	def get_beta_col(self,cat, ncat):
		print "Computing ellipticity-position misalignment angle between each object and its nearest neighbour."
		
		dx = cat["x_image"] - ncat["x_image"]
		dy = cat["y_image"] - ncat["y_image"]

		# position angle of the separation vector, in sky coordinates
		# has bounds [-pi,pi]
		theta = np.arctan(dy/dx)

		# position angle of the central galaxy in stamp coordinates
		eres = self.res["e1"]+ 1j*self.res["e2"]
		phi = np.angle(eres)/2



		# cos(beta) = Rneigh.ecent 
		# where Rneigh, ecent are unit vectors
		beta = (phi - theta)
		np.putmask(beta,np.invert(np.isfinite(beta)),0)

		# Impose bounds as above
		sel1=beta>np.pi/2
		sel2=beta<-1.0*np.pi/2
		beta[sel1] = np.pi/2 - (beta[sel1] - np.pi/2)
		beta[sel2] = -1.*np.pi/2 - (beta[sel2] + np.pi/2)

		beta/=np.pi

		return beta

		#self.res = arr.add_col(self.res,"dbeta",beta)

	def angle_cols(self, ncat):
		self.get_phi_col(ncat)
		self.get_beta_col(ncat)

	def get_2pt(self, corr1, corr2, nbins=12, error_type="bootstrap", xmin=1, xmax=300, units="arcmin"):

                import treecorr as tc
		correlation_lookup = {"psf":("mean_psf_e%d_sky", self.res), "gal":("e%d", self.res)}

		if hasattr(self,"truth"): correlation_lookup["int"] = ("intrinsic_e%d", self.truth)

		c1, data1 = correlation_lookup[corr1]
		c2, data2 = correlation_lookup[corr2]
		print "Will correlate columns %s and %s."%(c1,c2), 

		cat1 = tc.Catalog(ra=self.res["ra"]*60, dec=self.res["dec"]*60, ra_units=units, dec_units=units, g1=data1[c1%1], g2=data1[c1%2])
		cat2 = tc.Catalog(ra=self.res["ra"]*60, dec=self.res["dec"]*60, ra_units=units, dec_units=units, g1=data2[c2%1], g2=data2[c2%2])

		gg = tc.GGCorrelation(nbins=nbins, min_sep=xmin, max_sep=xmax, sep_units=units)

		gg.process(cat1,cat2)

		setattr(gg, "theta", np.exp(gg.logr))

		print "stored"
		self.corr[(corr1,corr2)]=gg

	def get_alpha_from_2pt(self):
		xi_pp = self.corr["psf","psf"].xip
		xi_gp = self.corr["gal","psf"].xip

		egal = self.res["e1"] + 1j*self.res["e2"]
		epsf = self.res["mean_psf_e1_sky"] + 1j*self.res["mean_psf_e2_sky"]

		alpha = ( xi_gp - np.mean(egal).conj()*np.mean(epsf) ) / ( xi_pp - abs(np.mean(epsf))*abs(np.mean(epsf)) )
		x = self.corr["psf","psf"].theta

		return x, alpha

	def get_fofz(self, nbins, error_type="bootstrap", zmax=1.8, binning="equal_number", T_B=False, return_number=False, split_type="im3shape", bin_type="mean_z_bpz"):
		"""Calculate the bulge fraction in bins using the given definition"""
		fr=[]
		e1=[]
		Tofz=[]
		bTofz=[]
		nofz_bulge=[]
		nofz_disc=[]

		if binning is "equal_number":
			bins = di.find_bin_edges(self.res[bin_type][self.res[bin_type]<zmax] , nbins)
		elif binning is "uniform":
			bins = np.linspace(0, zmax, nbins+1)

		for i, lower in enumerate(bins[:-1]):
			print i, 
			upper = bins[i+1]
			sel = (self.res[bin_type]>lower) & (self.res[bin_type]<upper)
			bulge_frac, selb = self.get_bulge_fraction(self.res[sel], split_type=split_type, return_mask=True)
			fr.append(bulge_frac)
			nofz_bulge.append(self.res[sel][selb].size)
			nofz_disc.append(self.res[sel][np.invert(selb)].size)

			if error_type is "bootstrap":
				e1.append(di.bootstrap_error(50, self.res[sel], self.get_bulge_fraction, additional_args=["split_type", "return_mask"], additional_argvals=[split_type, False]))
			else:
				e1.append(1.0 / np.sqrt(self.res[(sel)].shape[0]) )

			if T_B:
				Tofz.append(np.mean(self.res[(sel)]["T_B"]))
				bTofz.append(np.mean(self.res[(sel)]["T_B"]))

		z = (bins[1:]+bins[:-1])/2
		if T_B:
			out = (z, np.array(fr), np.array(e1), np.array(Tofz), np.array(bTofz))
		if return_number:
			out=(z,np.array(fr), np.array(e1), [nofz_disc, nofz_bulge])
		else:
			out = (z, np.array(fr), np.array(e1))

		return out


info_cuts =["(self.res['fails_unmasked_flux_frac']==0)", "(self.res['snr']>10)", "(self.res['snr']<10000)", "(self.res['mean_rgpp_rp']>1.1)", "(self.res['mean_rgpp_rp']<3.5)", "(self.res['radius']<5)", "(self.res['radius']>0.1)", "((self.res['ra_as']**2+self.res['dec_as']**2)**0.5<1.0)", "(self.res['chi2_pixel']>0.5)", "(self.res['chi2_pixel']<1.5)", "(self.res['min_residuals']>-0.2)", "(self.res['max_residuals']<0.2)", "(self.res['mean_psf_fwhm']<7.0)", "(self.res['mean_psf_fwhm']>0.0)", "(self.res['error_flag']==0)"]

class blank_stamp:
	def __init__(self, boxsize):
		self.array = np.zeros((boxsize,boxsize))

class dummy_im3shape_options():
	def __init__(self):
		self.psfex_rerun_version="y1a1-v02"
		self.verbosity=2
		rescale_stamp = "Y"
		self.upsampling = 5
		self.n_central_pixel_upsampling = 4
		self.n_central_pixels_to_upsample = 5
		self.padding = 4
		self.stamp_size = 48
		self.psf_input = "bundled-psfex"
		self.use_image_flags="Y"


class meds_wrapper(i3meds.I3MEDS):
	def __init__(self, filename, options=None,update=False, model="disc"):
		if options is None:
			options = p3s.Options("/global/cscratch1/sd/sws/hoopoe-image-simulations/end-to-end/end-to-end_code/%s_sim.ini"%model)
		super(meds_wrapper, self).__init__(filename, options)
		setattr(self, "filename", self._filename)

		if update:
			self._fits.close()
			print "Beware: FITS file can be overwritten in update mode."
			self._fits = fitsio.FITS(filename, "rw")

	def setup_dummy_im3shape_options(self):
		self.options = dummy_im3shape_options()
		
	def _get_extension_name(self, type):
		"""
		    Extended version of the function in meds.py.
		    Now can be used to get simulation specific cutouts.
		"""
		if type=='image':
			return "image_cutouts"
		elif type=="weight":
			return "weight_cutouts"
		elif type=="seg":
			return "seg_cutouts"
		elif type=="bmask":
			return "bmask_cutouts"
		elif type=="model":
			return "model_cutouts"
		elif type=="noise":
			return "noise_cutouts"
		elif type in [hdr.read_header()["EXTNAME"] for hdr in self._fits[1:]]:
			return type
		else: raise ValueError("bad cutout type '%s'" % type)

	def get_zpmag_scaling(self, iobj, iexp):

		exposure_id = self._cat["file_id"][iobj,iexp]
		zpmag = self._image_info["magzp"][exposure_id]
		zpmag_coadd = self._image_info["magzp"][0]

		return 10 **((zpmag_coadd - zpmag)/2.5)

	def get_stack_correction(self, iobj, silent=False):
		boxsize = self._cat["box_size"][iobj]
		images = np.array([f for f in self._cat["file_id"][iobj] if f>0])
		nexp = len(images)

		if not silent:
			print "Object %d has %d exposures"%(iobj, nexp)

		correction_stack = np.ones((boxsize*(nexp+1), boxsize))

		for iexp in xrange(nexp+1) :
			alpha = self.get_zpmag_scaling(iobj, iexp)
			correction_stack[boxsize*iexp : (boxsize*iexp + boxsize)] *= alpha

		return correction_stack

	def get_coadd_cat(self, path=None):
		if path is None:
			path = os.path.dirname(self._filename).split("meds")[0]
			
		tile = os.path.basename(self._filename)[:12]


	def remove_noise(self, silent=False, outdir="noisefree"):
		p0 = self._fits["image_cutouts"].read()
		real = self._fits["model_cutouts"]
		noise = self._fits["noise_cutouts"]

		if not silent:
			print "will remove noise"

		for iobj, object_id in enumerate(self._cat["id"]):
			if not silent:
				print iobj, object_id
			# Find the relevant index range for this object
			i0 = self._cat["start_row"][iobj][0]
			i1 = self._cat["start_row"][iobj][0]+self._cat["box_size"][iobj]*self._cat["box_size"][iobj]*self._cat["ncutout"][iobj]

			# COSMOS profile + neighbours
			scale_stack = self.get_stack_correction(iobj, silent=True).flatten()
			pixels0 = p0[i0:i1] - noise[i0:i1]*scale_stack
			pixels0*=p0[i0:i1].astype(bool).astype(int)

			import pdb ; pdb.set_trace()

			p0[i0:i1] = pixels


		print "Writing to MEDS file"
		self.clone(data=p0, colname="image_cutouts", newdir=outdir)

		return 0


	def remove_neighbours(self, silent=False, outdir="neighbourfree", noise=True, remove_masks=True):
		p0 = self._fits["image_cutouts"].read()
		real = self._fits["model_cutouts"]
		try:
			pn = self._fits["noise_cutouts"]
		except:
		    print "WARNING: No noise cutouts found"
		    pn = np.zeros_like(p0) 
		seg = self._fits["seg_cutouts"].read()

		if not silent:
			print "will remove neighbours"

		for iobj, object_id in enumerate(self._cat["id"]):
			if not silent:
				print iobj, object_id
			# Find the relevant index range for this object
			i0 = self._cat["start_row"][iobj][0]
			i1 = self._cat["start_row"][iobj][0]+self._cat["box_size"][iobj]*self._cat["box_size"][iobj]*self._cat["ncutout"][iobj]

			mask = seg[i0:i1]
			if remove_masks and (np.unique(mask).size>2):
				# Get the central pixel in the coadd seg map
				nx = self._cat["box_size"][iobj]
				#seg_id = mask[:nx*nx].reshape(nx,nx)[nx/2,nx/2]
				seg_id = self._cat["number"][iobj]
				# Remove any masks other than the one around the central object
				np.putmask(mask, mask!=seg_id, 0)

			scale_stack = self.get_stack_correction(iobj, silent=True).flatten()

			# COSMOS profile + noise
			if noise:
				neighbours = p0[i0:i1]- real[i0:i1]*scale_stack - pn[i0:i1]
			else:
				neighbours = p0[i0:i1]- real[i0:i1]*scale_stack
			pixels = p0[i0:i1] - neighbours
			pixels *= p0[i0:i1].astype(bool).astype(int)

			p0[i0:i1] = pixels
			seg[i0:i1] = mask

		print "Writing to MEDS file"
		self.clone(data=p0, colname="image_cutouts", data2=seg, colname2="seg_cutouts", newdir=outdir, rename=["seg_cutouts", "hoopoe_seg_cutouts"])

		return 0

	def remove_model_bias(self, shapecat, silent=False, outdir="analytic", noise=True, neighbours=True, remove_masks=False):
		p0 = self._fits["image_cutouts"].read()
		real = self._fits["model_cutouts"]
		seg = self._fits["seg_cutouts"].read()

		if not silent:
			print "will insert analytic profiles"

		if not hasattr(shapecat, "res"):
			setattr(shapecat, "res", shapecat.truth)

		for iobj, object_id in enumerate(self._cat["id"]):
			if not silent:
				print iobj, object_id
			# Find the relevant index range for this object
			i0 = self._cat["start_row"][iobj][0]
			i1 = self._cat["start_row"][iobj][0]+self._cat["box_size"][iobj]*self._cat["box_size"][iobj]*self._cat["ncutout"][iobj]

			# First subtract off the input profile flux
			existing = real[i0:i1]
			
			pixels = p0[i0:i1]-existing
			pixels*=p0[i0:i1].astype(bool).astype(int)

			icat = np.argwhere(shapecat.res["id"]==object_id)
			if len(icat)>0:
				icat=icat[0,0]
				res = shapecat.res[icat]
			else:
				print "Warning: no im3shape results found for object."
				res = np.random.choice(shapecat.res)
				icat = np.argwhere(shapecat.res["id"]==res["id"])[0,0]

			if hasattr(shapecat, "truth"):
				truth= shapecat.truth[icat]

			gal = self.construct_analytic_profile(res, shapecat.fit, truth=truth)
			stack = self.make_stack(gal,iobj)

			mask = seg[i0:i1]
			if remove_masks and (np.unique(mask).size>2):
				# Get the central pixel in the coadd seg map
				nx = self._cat["box_size"][iobj]
				seg_id = mask[:nx*nx].reshape(nx,nx)[nx/2,nx/2]
				# Remove any masks other than the one around the central object
				np.putmask(mask, mask!=seg_id, 0)

			scale_stack = self.get_stack_correction(iobj, silent=True).flatten()

			if noise and (not neighbours):
				remove = p0[i0:i1]- real[i0:i1] * scale_stack - noise[i0:i1] * scale_stack
			elif (not noise) and (not neighbours):
				remove = p0[i0:i1]- real[i0:i1] * scale_stack
			elif (not noise) and neighbours:
				remove = noise * scale_stack
			elif noise and neighbours:
				remove = 0.0

			# noise + neighbours + analytic profile
			try:
				pixels = pixels + stack - remove
			except:
				import pdb ; pdb.set_trace()

			p0[i0:i1] = stack

			

		print "Writing to MEDS file %s"%outdir
		self.clone(data=p0, colname="image_cutouts", newdir=outdir)

		return 0

	def clone(self, data=None, colname="image_cutouts", data2=None, colname2="image_cutouts", newdir=None, rename=[None,None]):
		out=fitsio.FITS("%s/%s"%(newdir,os.path.basename(self._filename).replace(".fz","")), "rw")

		if rename[0] is not None:
			print "Will rename column %s as %s"%(rename[0], rename[1])
			dat = self._fits[rename[0]].read()
			hdr = self._fits[rename[0]].read_header()

			out.write(dat, header=hdr)
			out[-1].write_key('EXTNAME', rename[1])

		for j, h in enumerate(self._fits[1:]):
			hdr=h.read_header()
			extname=hdr["extname"]
			if extname in ["model_cutouts", "noise_cutouts"]:
				continue
			print extname
			if (data is not None) and extname==colname:
				print "new data in HDU %s"%colname
				dat = data
			elif (data2 is not None) and extname==colname2:
				print "new data in HDU %s"%colname2
				dat = data2
			else:
				dat=h.read()
			out.write(dat, header=hdr)
			out[-1].write_key('EXTNAME', extname)

		out.close()

		print "Compressing output file"
		filename = "%s/%s"%(newdir,os.path.basename(self._filename).replace(".fz",""))
		fpack_cmd = "fpack  -qz 4.0 -t 10240,1 {}".format(filename)
		os.system(fpack_cmd)
		os.system("rm -rf %s"%filename)

	def make_stack(self,profile,iobj):
		boxsize = self._cat["box_size"][iobj]
		nexp = self._cat["ncutout"][iobj]

		stack = [np.zeros(boxsize*boxsize)]

		for iexp in xrange(1,nexp):
			wcs, offset = self.get_wcs(iobj,iexp)

			# Get the relevant zero point magnitude scaling
			i_image = self._cat["file_id"][iobj][iexp]
			magzp = self._image_info[i_image]["magzp"]
			rf = 10**( (magzp - 30) / 2.5 )

			psf = self.get_bundled_psfex_psf(iobj, iexp, return_profile=True)

			if (profile is not None) and (psf is not None):
				final = galsim.Convolve([profile,psf[0]]) * rf
			elif profile is None:
				final = psf[0]
			elif psf is None:
				final = profile
			
			if final is not None:
				stamp = galsim.ImageD(boxsize, boxsize)
				tmp = final.drawImage(stamp, scale=1, offset=offset, method='no_pixel')
			else:
				stamp=blank_stamp(boxsize)

			stack.append(stamp.array.flatten())

		return np.concatenate(stack)


#	def get_bundled_psfex_psf(self, iobj, iexp, return_profile=False):
#
#		wcs_path = self.get_source_info(iobj,iexp)[3].strip()
#
#		wcs_path = check_wcs(wcs_path)
#
#		#Find the exposure name
#		source_path = self.get_source_path(iobj, iexp)
#		source = os.path.split(source_path)[1].split('.')[0]
#
#		#PSF blacklist
#		if source in self.blacklist:
#			if self.options.verbosity>2:
#				print "%s blacklisted" % source
#			return None
#
#		#get the HDU name corresponding to that name.
#		#Might need to tweak that naming
#		hdu_name = "psf_"+self.options.psfex_rerun_version+"_"+source+"_psfcat"
#
#		#Maybe we have it in the cache already
#		if hdu_name in PSFEX_CACHE:
#			psfex_i = PSFEX_CACHE[hdu_name]
#			if psfex_i is None:
#				return None
#
#		else:
#    		#Turn the HDU into a PSFEx object
#    		#PSFEx expects a pyfits HDU not fitsio.
#    		#This is insane.  I know this.
#			import galsim.des
#			try:
#				pyfits = galsim.pyfits
#			except AttributeError:
#				from galsim._pyfits import pyfits
#
#			try:
#				hdu = pyfits.open(self._filename)[hdu_name]
#			except KeyError:
#				PSFEX_CACHE[hdu_name] = None
#				if self.options.verbosity>3:
#					print "Could not find HDU %s in %s" % (hdu_name, self._filename)
#
#				return None
#
#			try:
#				psfex_i = galsim.des.DES_PSFEx(hdu, wcs_path)
#			except IOError:
#				print "PSF bad but not blacklisted: %s in %s"%(hdu_name, self._filename)
#				psfex_i = None
#
#			PSFEX_CACHE[hdu_name] = psfex_i
#
#		if psfex_i == None:
#			return None
#
#		#Get the image array
#		return self.extract_psfex_psf(psfex_i, iobj, iexp, return_profile=return_profile, orig_func=True)
#
#
#
#	def extract_psfex_psf(self, psfex_i, iobj, iexp, return_profile=False, orig_func=False):
#		psf_size=(self.options.stamp_size+self.options.padding)*self.options.upsampling
#		orig_col = self['orig_col'][iobj][iexp]
#		orig_row = self['orig_row'][iobj][iexp]
#
#		x_image_galsim = orig_col+1
#		y_image_galsim = orig_row+1
#
#		if orig_func:
#			psf = utils.getPSFExarray(psfex_i, orig_col, orig_row, psf_size, psf_size, self.options.upsampling, return_profile=return_profile)
#			return psf
#
#		else:
#			psf = psfex_i.getPSF(galsim.PositionD(x_image_galsim, y_image_galsim))
#
#			image = galsim.ImageD(psf_size, psf_size)
#
#			if return_profile:
#				return psf

	def get_wcs(self, iobj, iexp):

		wcs_path = self.get_source_path(iobj, iexp)

		wcs_path = check_wcs(wcs_path)

		orig_col = self['orig_col'][iobj][iexp]
		orig_row = self['orig_row'][iobj][iexp]        
		image_pos = galsim.PositionD(orig_col,orig_row)
		wcs = galsim.FitsWCS(wcs_path)

		ix = int(math.floor(image_pos.x ))
		iy = int(math.floor(image_pos.y ))

		# Offset of the galaxy centroid from the stamp centre
		offset = (image_pos - galsim.PositionD(ix,iy))
		offset= galsim.PositionD(offset.x+0.5, offset.y+0.5)

		return wcs.local(image_pos), offset


	def construct_analytic_profile(self, res, fittype, truth=None):
		if fittype in sersic_indices.keys():
			n = sersic_indices[fittype]
		elif fittype=="bord":
			n = 1 + res["is_bulge"].astype(int)*3
			fittype = {4: "bulge", 1: "disc"} [n] 

		if truth is not None:
			e1 = truth["intrinsic_e1"]
			e2 = truth["intrinsic_e2"]
			g1 = truth["true_g1"]
			g2 = truth["true_g2"]
			hlr = truth["hlr"]
			mag = truth["mag"]
			flux = 10 **(-1.0 * (mag - 30) / 2.5)

			if (g1==-9999.) or (g2==-9999.):
				return None
		else:
			e1 = res["e1"][0,0]
			e2 = res["e2"][0,0]
			g1 = 0
			g2 = 0
			flux = res["%s_flux"%fittype][0,0]

		gal=galsim.Sersic(n=n,half_light_radius=hlr)

		gal=gal.withFlux(flux)

		e1 += g1
		e2 += g2
		if e1>1: e1=1
		if e1<-1: e1=-1
		if e2>1: e2=1
		if e2<-1: e2=-1

		try:
			shear = galsim.Shear(g1=e1, g2=e2)
			print "applying shear g1=%1.3f,g2=%1.3f"%(e1,e2)
		except:
			print "ERROR: unphysical shear"
			shear= galsim.Shear(g1=0, g2=0)

		gal.shear(shear)

		return gal

	def download_source_data(self):

		import pdb ; pdb.set_trace()
		
		for iobj, object_id in enumerate(self._cat["id"]):
			nexp = self._cat["ncutout"][iobj]
			i0 = 0
			i1 = nexp
			if not coadd:
				i0+=1
			if not red:
				i1=0

	def download_source_images(self, coadd=False, red=True):
		
		for iobj, object_id in enumerate(self._cat["id"]):
			nexp = self._cat["ncutout"][iobj]
			i0 = 0
			i1 = nexp
			if not coadd:
				i0+=1
			if not red:
				i1=0

			for iexp in xrange(i0,i1):
				file_id = self._cat["file_id"][iobj, iexp]
				path = self._image_info["image_path"][file_id].strip()
				path = check_wcs(path,iexp)

	def i3s(self, iobj, sub_image=[None,None,None], weights=[None,None,None], show=False, save=None, ncols=1, col=1, coadd_objects_id=-9999, return_vals=False):
		# Setup im3shape inputs
		opt = self.options
		if coadd_objects_id!=-9999:
			iobj=np.argwhere(self._fits["object_data"].read()["id"]==coadd_objects_id)[0,0]
			print coadd_objects_id, iobj
		inputs = self.get_im3shape_inputs(iobj)
		psfs = inputs.all('psf')
		bands = self.get_band() * (len(inputs))

		galaxy_stack = inputs.all('image')

		if weights[0] is None:
			weights_stack = inputs.all('weight')
		else:
			weights_stack = weights

		# Run the MCMC
		if sub_image[0] is not None:
			print "Warning: will substitute in new image (weight, PSFs, masks upchanged)"
			galaxy_stack = sub_image
			result, best_img, images, weights = p3s.analyze_multiexposure( sub_image, psfs, weights_stack, inputs.all('transform'), self.options, ID=3000000, bands=bands)
		else:
			print "bands:", bands
			result, best_img, images, weights = p3s.analyze_multiexposure( inputs.all('image'), psfs, weights_stack, inputs.all('transform'), self.options, ID=coadd_objects_id, bands=bands)

		unmasked_fraction = [ p3s.utils.get_unmasked_flux_fraction(image,wt) for i,(image, wt) in enumerate(zip(images, weights_stack)) ]
		print "Mean Mask Fraction : %f"%(1-np.array(unmasked_fraction).mean())
		boxsize = wt.shape[0]
		nexp = len(weights_stack)
		effective_npix = boxsize*boxsize * nexp * np.array(unmasked_fraction).mean()
		print "Effective npix : boxsize^2 x Num Exp x Unmasked Pixel Frac = %f^2 x %f x %f = %f"%(boxsize, nexp, np.array(unmasked_fraction).mean(), effective_npix)
		print "chi2 per pixel : -2 x %f / %f = %f"%(result.obj.likelihood, effective_npix, -2*result.obj.likelihood / effective_npix)

		if isinstance(save,str):
			galaxy_stack = inputs.all('image')
			wt = inputs.all('weight')
			image = np.hstack(tuple(galaxy_stack))
			np.savetxt(save+"/model.txt", best_img )
			np.savetxt(save+"/image.txt", image )
			np.savetxt(save+"/weights.txt", np.hstack(tuple(wt)) )

		image = np.hstack(tuple(galaxy_stack)) 
		wt = inputs.all('weight')

	
		if show:
			import pylab as plt
			cmap = plt.get_cmap("jet")

			nrow = 2

			plt.subplot(nrow,ncols, col)
			plt.title("Image")
			plt.imshow(image, interpolation="none", cmap=cmap)
			plt.colorbar()
			plt.subplot(nrow,ncols,ncols+col)
			plt.title("Best fit")
			plt.imshow(best_img, interpolation="none", cmap=cmap)
			print "Pixel max:%f"%best_img.max()
			plt.colorbar() 

	
		if return_vals:
			return result.get_params(), image, best_img, np.hstack(tuple(wt)), inputs.all('transform')
		else:
			return result, result.get_params()



def check_wcs(wcs_path, iexp=-9999):
	if not os.path.exists(wcs_path):
		wcs_path = "/share/des/disc3/samuroff/y1/data/OPS"+wcs_path.split("OPS")[-1].strip()
		if iexp==0:
			wcs_path=wcs_path.replace(".fz","")

	if not os.path.exists(wcs_path):
		print "Warning. Could not find WCS, so will download it to %s."%wcs_path
		print "Are you sure this is what you want to do?"
		source = wcs_path.replace("/share/des/disc3/samuroff/y1/data/OPS/", "/global/cscratch1/sd/sws/data/OPS/")
		os.system("mkdir -p %s"%os.path.dirname(wcs_path))
		os.system("scp sws@edison.nersc.gov:%s %s"%(source,os.path.dirname(wcs_path)))

	return wcs_path




#& (self.res['chi2_pixel']>0.5) &



def add_col(rec, name, arr=[], dtype=None):
	"""Generic function to add a new column to a structured numpy array.
	Borrows heavily from Tomek's code."""

	if len(arr)==0:
		arr=np.zeros(len(rec))

	arr = np.asarray(arr)
	if dtype is None:
		dtype = arr.dtype

	newdtype = np.dtype(rec.dtype.descr + [(name, dtype)])
	newrec = np.empty(rec.shape, dtype=newdtype)
	for field in rec.dtype.fields:
		newrec[field] = rec[field]

	newrec[name] = arr

	return newrec

def split_by(array, column_name, pivot, return_mask=False, logger=None):
	"""Return halves of a structured array split about a given value in one of the columns."""

	lower= array[(array[column_name]<pivot)]
	upper= array[(array[column_name]>pivot)]

	f = 1.0*len(upper)/len(array)

	if logger is not None:
		logger.info("Splitting array about %s=%3.4f (%3.2f percent above)"%(column_name,pivot, f*100.))

	if return_mask:
		return f, (array[column_name]<pivot)

	else:
		return f, upper, lower

def translate_infocut(info_val):
	bits_set = [i for i,j in enumerate(' '.join(word[::-1] for word in bin(info_val)[2:].split()) ) if int(j)==1]

	for bit_set in bits_set: 
		print bit_set, "%s : %s"%INFO_FLAGS[bit_set]

INFO_FLAGS={ 
0   : ("INFO_GOLD_MASK",  "This area is masked out in the gold catalog" ),
1   : ("INFO_BAD_MASK",  "This object is flagged in the gold catalog" ),
2   : ("INFO_MODEST",  "This modest classifier suggests this is not a galaxy" ),
3   : ("INFO_MASK",  "mask fraction > 0.75" ),
4   : ("INFO_EVALS",  "levmar_like_evals > 10000" ),
5   : ("INFO_SEXR_NEIGHBOURS",  "r-band sextractor flagged with 0x1, indicating bright neigbours" ),
6   : ("INFO_SEXR_BLEND",  "r-band sextractor flagged with 0x2, indicating blending" ),
7   : ("INFO_MASKFLUX",  "More than 25% of flux masked" ),
8   : ("INFO_SMALL_SNR",  "S/N < 10" ),
9   : ("INFO_HUGE_SNR",  "snr > 10000" ),
10   : ("INFO_SMALL_RGPP_RP",  "rgpp_rp < 1.1" ),
11   : ("INFO_LARGE_RGPP_RP",  "rgpp_rp > 3.5 (very large galaxy)" ),
12   : ("INFO_LARGE_RADIUS",  "radius > 5 arcsec" ),
13   : ("INFO_SMALL_RADIUS",  "radius < 0.1 arcsec" ),
14   : ("INFO_LARGE_CENTROID_SHIFT",  "centroid more than 1 arcsec from nominal" ),
15   : ("INFO_SMALL_CHI2",  "Chi2 per effective pixel < 0.5" ),
16   : ("INFO_LARGE_CHI2",  "Chi2 per effective pixel > 1.5" ),
17   : ("INFO_MINRES",  "Normed residuals < -0.2 somewhere" ),
18   : ("INFO_MAXRES",  "Normed residuals > 0.2 somewhere" ),
19   : ("INFO_LARGE_PSF",  "Very large PSF" ),
20   : ("INFO_NEG_PSF",  "Negative PSF FWHM" ),
21   : ("INFO_ERROR",  "One or more error flags is set" ) }

