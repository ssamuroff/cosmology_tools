import numpy as np
import pyfits as pf
import fitsio as fio
import glob
import tools.diagnostics as di
import pylab as plt
import tools.arrays as arr
from plots import im3shape_results_plots as i3s_plots

names = {"snr": "snr", "rgp" : "mean_rgpp_rp", "e1": "e1", "e2" : "e2", "iterations":"levmar_iterations", "psf_e1" : "mean_psf_e1", "psf_e2" : "mean_psf_e2", "psf_fwhm": "mean_psf_fwhm", "nflag": "neighbour_flag"}

class shapecat(i3s_plots):
	"""Class for handling im3shape results and,
	   if relevant, truth catalogues.
	"""
	def __init__(self, res=None, truth=None, coadd=False):
		if res is None and truth is None:
			print "Please specify at least one of res=, truth="
		self.res_path = res
		self.truth_path = truth
		self.coadd_path= coadd

		print "Initialised."

	def load_from_array(self, array, name="res"):
		setattr(self, name, array)

	def load(self, res=True, truth=False, coadd=False, postprocessed=True, keyword="DES", apply_infocuts=True, ext=".fits"):
		
		if "%s"%ext in self.res_path:
			files=[self.res_path]
		else:
			files = glob.glob("%s/*.%s"%(ext,self.res_path))
		single_file=False
		print "loading %d results file(s) from %s"%(len(files),self.res_path)
		if len(files)>1:
			self.res, self.files, i = di.load_results(res_path =self.res_path, format=ext, apply_infocuts=apply_infocuts, keyword=keyword, postprocessed=postprocessed, return_filelist=True)
		else:
			if ext.lower()==".fits":
				self.res = fio.FITS(files[0])[1].read()
			elif ext.lower()==".txt":
				self.res = np.genfromtxt(files[0], names=True)
			self.files=files
			single_file=True
			i=None


		if truth:
			if ".fits" in self.truth_path:
				files=[self.truth_path]
			else:
				files = glob.glob("%s/*.fits"%self.truth_path)
			single_file=False

			print "loading truth files from %s"%self.truth_path
			if len(files)>1:
				self.truth = di.load_truth(truth_path=self.truth_path, match=self.files, res=self.res)
			else:
				self.truth = pf.getdata(files[0])

		if truth and res:
			self.res, self.truth = di.match_results(self.res,self.truth)
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

	def add_bpz_cols(self, fil="/home/samuroff/NSEVILLA_PHOTOZ_TPL_Y1G103_1_bpz_highzopt_2_9_16.fits", array=None, exclude=None):
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

	def get_positions(self, postype):
		"""Extract the position coordinates of objects, either from a truth
		   table or from an im3shape results table."""
		if postype not in ["world", "pixel"]:
			"Please choose a coordinate system ('world' or 'pixel')"

		if postype=="world":
			xname,yname = "ra","dec"
		elif postype=="pixel":
			xname,yname = "X_IMAGE","Y_IMAGE"

		print "Obtained position data from",

		try:
			x,y = self.truth[xname],self.truth[yname]
			print "truth table"
		except:
			try:
				x,y = self.res[xname],self.res[yname]
				print "results table"
			except:
				x,y = self.coadd[xname],self.coadd[yname]
				print "coadd catalogue"

		return x,y


	def hist(self, param, normed=1, histtype="step", alpha=1.0, truth=False, res=True, nbins=25, xlim_upper=None, xlim_lower=None, colour="purple", infoflags=False, errorflags=False, neighbourflags=False, starflags=False, linestyle="solid", extra_cuts=None):
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
		plt.hist(data[colname][sel], bins=bins, histtype=histtype, color=colour, normed=normed, lw=2.0, linestyle=linestyle, alpha=alpha)

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

	def get_infocuts(self, exclude="None"):
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
		return cuts

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

	def get_neighbours(self, truth=True):
		import copy
		ni=np.zeros_like(self.res["e1"])
		for i,j in enumerate(self.truth["nearest_neighbour_index"]):
			if j in self.res["id"]:
				ni[i]=np.argwhere(self.res["id"]==j)[0,0]

		neighbour_cat = copy.deepcopy(self)
		setattr(neighbour_cat,"res",self.res[ni.astype(int)])
		setattr(neighbour_cat,"truth",self.truth[ni.astype(int)])

		return neighbour_cat

	def get_phi_col(self,neighbour_cat):
		print "Computing misalignment angle between each object and its nearest neighbour."
		de1 = self.res["e1"] - neighbour_cat.res["e1"]
		de2 = self.res["e2"] - neighbour_cat.res["e2"]
		dphi = 0.5 * np.arctan2(de2, de1)
		self.res = arr.add_col(self.res,"dphi",dphi)

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





		