import numpy as np
import scipy as sp
import astropy.io.fits as pyfits
import os, pdb, copy
import tools.diagnostics as di
import tools.arrays as arr
import tools.plots as plots
import tools.shapes as sh
import tktools
import gc
from scipy.interpolate import Rbf

import scipy.optimize as opt
import fitsio as fi

import matplotlib.colors
import matplotlib
import pylab as plt
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['font.size']=14
matplotlib.rcParams['legend.fontsize']=14
matplotlib.rcParams['xtick.major.size'] = 10.0
matplotlib.rcParams['ytick.major.size'] = 10.0

class nbc(plots.im3shape_results_plots, sh.shapecat):
	def __init__(self, from_disc=False, from_memory=True, res_dirname=None, truth_dirname=None, res=None, truth=None):
		if from_memory:
			print "initialising from memory"
			self.res = res
			self.truth = truth
		elif from_disc:
			print "initialising from disc"
			self.results_path = res_dirname
			self.truth_path = truth_dirname
			self.load(truth=True)
		else:
			print "Please specify either from_disc=True or from_memory=True"

	def get_split_data(self, cat, weights=None, method="random"):
		self.res = cat.res
		# Fix this so the random split for a given input is always the same
		np.random.seed(9000)
		if method=="random":
			print "Subdividing catalogue by coadd objects ID"
			sel = np.random.randint(0,2,cat.res.size).astype(bool)
		elif method=="cosmos":
			print "Subdividing catalogue by COSMOS identifier"
			cosmos_ids = np.unique(cat.truth["cosmos_ident"])
			selected_ids = np.random.choice(cosmos_ids, cosmos_ids.size/2, replace=False)
			sel = np.in1d(cat.truth["cosmos_ident"], selected_ids)

		self.res1 = self.res[sel]
		self.res2 = self.res[np.invert(sel)]

		if hasattr(cat,"truth"):
			self.truth = cat.truth
			self.truth1 = cat.truth[sel]
			self.truth2 = cat.truth[np.invert(sel)]

		if weights is not None:
			return weights[sel], weights[np.invert(sel)]

	def bias_by_grid_cells(self, bias_name, catalogue_to_calibrate):
		print "Will calibrate in grid cells"
	
		print "Found grid points of dimensions", self.bias_grid.shape

		bias = np.zeros(catalogue_to_calibrate.size)

		# Go through each size bin in turn
		for i, rbounds in enumerate(zip( np.unique(self.bias_grid["rgp_lower"]), np.unique(self.bias_grid["rgp_upper"]) )):
			r_lower = rbounds[0]
			r_upper = rbounds[1]

			# Make the selection within these size bounds
			select_r = (catalogue_to_calibrate["mean_rgpp_rp"]<r_upper) & (catalogue_to_calibrate["mean_rgpp_rp"]>r_lower)

			# Get the SNR edges within this size bin
			snr_lower = self.bias_grid["snr_lower"][(self.bias_grid["i"]==i)]
			snr_upper = self.bias_grid["snr_upper"][(self.bias_grid["i"]==i)]

			sub_array = copy.copy(catalogue_to_calibrate[select_r])
			sub_bias = copy.copy(bias[select_r])

			for j, sbounds in enumerate(zip(snr_lower, snr_upper)):
				s_lower = sbounds[0]
				s_upper = sbounds[1]

				print "%d [%3.2f-%3.2f] %d [%3.2f-%3.2f]"%(i, s_lower, s_upper, j, r_lower, r_upper ),

				# Make the SNR selection
				select_s = (sub_array["snr"]<s_upper) & (sub_array["snr"]>s_lower)

				if bias_name=="m":
					# Single value in this grid cell
					val = self.bias_grid["m"][(self.bias_grid["i"]==i) & (self.bias_grid["j"]==j)][0]
				elif bias_name=="a":
					# Single value in this grid cell
					val = self.bias_grid["alpha"][(self.bias_grid["i"]==i) & (self.bias_grid["j"]==j)][0]

				sub_array[select_s] = val 

				print " %d galaxies, %s=%2.3f"%(sub_bias[select_s].size, bias_name, val)

			bias[select_r] = sub_array


		return bias
					

	def apply(self, names=["m","a"], split_half=2, scheme="rbf"):

		if scheme=="rbf":
			print "Using RBF interpolation"
			use_rbf=True
		elif scheme=="grid":
			print "Using simple gridded bias"
		else:
			print "Using polynomial fit"

		# Choose input data
		if split_half>0:
			exec "catalogue_to_calibrate = self.res%d"%split_half
		else:
			catalogue_to_calibrate = self.res

		# Apply the best-fit calibration coefficients

		for bias_name in names:
			gc.collect()
			print "creating column %s"%bias_name
			if scheme=="polynomial":
				com2 = "a0 "+", a%d "*self.optimised_coefficients_m[1:].size +"=tuple(self.optimised_coefficients_%s)"%bias_name
				com2 = com2%tuple(np.linspace(1,17,17))
				exec com2

				com="eval_%s(np.array([catalogue_to_calibrate['snr'],catalogue_to_calibrate['mean_rgpp_rp']]).T"%bias_name + ", a%d "*self.optimised_coefficients_m.size+")"
				com=com%tuple(np.linspace(0,17,18))
				exec "bias=%s"%com
			elif scheme=="grid":
				bias = self.bias_by_grid_cells(bias_name, catalogue_to_calibrate)
			elif use_rbf:
				bias = self.do_rbf_interpolation(bias_name, catalogue_to_calibrate)
				

			if split_half>0:
				if split_half==1:
					self.res1 = arr.add_col(self.res1, bias_name, np.zeros_like(bias))
					self.res1[bias_name] = bias
					if bias_name=="a":
						self.res1 = arr.add_col(self.res1, 'c1', bias*self.res1['mean_hsm_psf_e1_sky'])
						self.res1 = arr.add_col(self.res1, 'c2', bias*self.res1['mean_hsm_psf_e2_sky'])
						self.res1['c1'] = bias * self.res1['mean_hsm_psf_e1_sky']
						self.res1['c2'] = bias * self.res1['mean_hsm_psf_e2_sky']
				if split_half==2:
					self.res2 = arr.add_col(self.res2, bias_name, np.zeros_like(bias))
					self.res2[bias_name] = bias
					if bias_name=="a":
						self.res2 = arr.add_col(self.res2, 'c1', bias*self.res2['mean_hsm_psf_e1_sky'])
						self.res2 = arr.add_col(self.res2, 'c2', bias*self.res2['mean_hsm_psf_e2_sky'])
						self.res2['c1'] = bias * self.res2['mean_hsm_psf_e1_sky']
						self.res2['c2'] = bias * self.res2['mean_hsm_psf_e2_sky']
				
			else:
				if bias_name not in self.res.dtype.names:
					self.res = arr.add_col(self.res, bias_name, bias)
				else:
					self.res[bias_name]=bias
				if bias_name=="a":
					if "c1" not in self.res.dtype.names:
						self.res = arr.add_col(self.res, "c1", bias*self.res["mean_hsm_psf_e1_sky"])
						self.res = arr.add_col(self.res, "c2", bias*self.res["mean_hsm_psf_e2_sky"])
					else:
						self.res["c1"] = bias * self.res["mean_hsm_psf_e1_sky"]
						self.res["c2"] = bias * self.res["mean_hsm_psf_e2_sky"]

			print "Finished calibration"

	def compute(self, split_half=0, redshift_bin=None, weights=None, fit="bord", reweight_per_bin=False, resample_per_bin=False, apply_calibration=False, refdata=None, table_name=None, ellipticity_name="e", sbins=10, rbins=5, binning="equal_number",rlim=(1,3), slim=(10,1000)):
		print 'measuring bias'



		if split_half>0:
			exec "data = self.res%d"%split_half
			print "using half %d of the catalogue (%d objects)"%(split_half,data.size)
			if hasattr(self, "truth"):
				exec "tr = self.truth%d"%split_half
		else:
			data = self.res
			print "using the full catalogue (%d objects)"%(data.size)
			if hasattr(self, "truth"):
				tr = self.truth

		if redshift_bin is not None:
			print "Applying additional cuts"
			if weights is not None:
				weights = weights[data["des_bin"]==redshift_bin]
			tr = tr[data["des_bin"]==redshift_bin]
			data = data[data["des_bin"]==redshift_bin]

			print "bin contains %d galaxies"%data.size

		if fit.lower()!="bord":
			print "Using %s only galaxies"%fit
			val = int(fit.lower()=="bulge")
			sel_fit = data["is_bulge"]==val
		else:
			sel_fit = np.ones_like(data).astype(bool)



		data = data[sel_fit]
		tr = tr[sel_fit]
		wts=weights[sel_fit]

		sel_lim = (data["snr"]>slim[0]) & (data["snr"]<slim[1]) & (data["mean_rgpp_rp"]>rlim[0]) & (data["mean_rgpp_rp"]<rlim[1])
		data = data[sel_lim]
		tr = tr[sel_lim]
		wts = wts[sel_lim]

		if refdata is not None:
			import tools.likelihoods as lk
			if fit.lower()!="bord":
				print "Using %s only galaxies"%fit
				val = int(fit.lower()=="bulge")
				sel_fit = refdata.res["is_bulge"]==val
			else:
				sel_fit = np.ones_like(refdata.res["e1"]).astype(bool)
			refdata.res = refdata.res[sel_fit]
			edat = np.sqrt(refdata.res["e1"]**2+refdata.res["e2"]**2)
			eh = np.sqrt(data["e1"]**2+data["e2"]**2)

		if isinstance(binning,str) : 
			if binning.lower()=="uniform":
				snr_edges = np.logspace(np.log10(slim[0]),np.log10(slim[1]),sbins+1)
				rgp_edges = np.linspace(rlim[0],rlim[1],rbins+1)
			elif binning.lower()=="equal_number":
				snr_edges = di.find_bin_edges(np.log10(data["snr"]), sbins)
				rgp_edges = di.find_bin_edges(np.log10(data["mean_rgpp_rp"]), rbins)

			override_bin_edges = False
			nrbin = rbins
			nsbin = sbins
		elif isinstance(binning, tuple):
			snr_edges = binning[1]
			rgp_edges_low = np.unique(np.array(binning[0]).T[0])
			rgp_edges = np.hstack(( np.unique(np.array(binning[0]).T[0]), np.atleast_1d(np.unique(np.array(binning[0]).T[1])[-1]) ))
			nrbin = np.unique(np.array(binning[0]).T[0]).size
			override_bin_edges = True

		list_bias = []
		bias_grid=[]

#		b = di.get_bias(tr, data, nbins=5, apply_calibration=apply_calibration, ellipticity_name=ellipticity_name, binning="equal_number", names=["m","c","m11","m22","c11","c22"], silent=True)
#		print "Global biases:"
#		print "m11 : ", b["m11"]
#		print "m22 : ", b["m22"]
#		print "m : ", b["m"]
#		print "c11 : ", b["c11"]
#		print "c22 : ", b["c22"]
#		print "c : ", b["c"]

		if not override_bin_edges:
			print "Will do dynamic binning in SNR"

		for i in xrange(nrbin):
			if not override_bin_edges:	
				snr_samp = data["snr"][(np.log10(data['mean_rgpp_rp']) > rgp_edges[i]) & (np.log10(data['mean_rgpp_rp']) < rgp_edges[i+1])]
				snr_edges=di.find_bin_edges(np.log10(snr_samp), sbins)
			else:
				select_edges = (np.array(binning[0]).T[0]==rgp_edges[i])
				snr_edges_low = np.array(binning[1]).T[0][select_edges] 
				snr_edges_high = np.array(binning[1]).T[1][select_edges]
				snr_edges = np.hstack((snr_edges_low, np.atleast_1d(snr_edges_high[-1])))
				sbins = snr_edges_high.size
			for j in xrange(sbins):
				empty=False
				

				# Select in bins of snr and size
				# This is horrible. I know it, and I'm sorry.
				if not override_bin_edges:
					slow,shigh = 10**snr_edges[j], 10**snr_edges[j+1]
					rlow,rhigh = 10**rgp_edges[i], 10**rgp_edges[i+1]
					select = (data['snr'] > 10**snr_edges[j]) & (data['snr'] < 10**snr_edges[j+1]) & (data['mean_rgpp_rp'] > 10**rgp_edges[i]) & (data['mean_rgpp_rp'] < 10**rgp_edges[i+1])
					print "bin %d %d  snr = [%2.3f-%2.3f] rgpp/rp = [%2.3f-%2.3f]"%(j, i, 10**snr_edges[j], 10**snr_edges[j+1], 10**rgp_edges[i], 10**rgp_edges[i+1] )
				else:
					slow,shigh = 10**snr_edges_low[j], 10**snr_edges_high[j]
					rlow,rhigh = 10**binning[0][i][0], 10**binning[0][i][1]
					select = (data['snr'] > 10**snr_edges_low[j]) & (data['snr'] < 10**snr_edges_high[j]) & (data['mean_rgpp_rp'] > 10**binning[0][i][0]) & (data['mean_rgpp_rp'] < 10**binning[0][i][1])
					print "bin %d %d  snr = [%2.3f-%2.3f] rgpp/rp = [%2.3f-%2.3f]"%(j, i, 10**snr_edges_low[j], 10**snr_edges_high[j], 10**binning[0][i][0], 10**binning[0][i][1] )
				ngal = np.nonzero(select.astype(int))[0].size

				if resample_per_bin:
					select_data = (refdata.res["snr"]>slow) & (refdata.res["snr"]<shigh) & (refdata.res["mean_rgpp_rp"]>rlow) & (refdata.res["mean_rgpp_rp"]<rhigh) 
					subsample = di.get_selection_to_match(edat[select_data],eh[select],nbins=13)
				else:
					subsample = np.ones_like(data["e1"][select]).astype(bool)

				if reweight_per_bin:
					select_data = (refdata.res["snr"]>slow) & (refdata.res["snr"]<shigh) & (refdata.res["mean_rgpp_rp"]>rlow) & (refdata.res["mean_rgpp_rp"]<rhigh) 
					bin_wts=di.get_weights_to_match(edat[select_data],eh[select],nbins=15)
				else:
					bin_wts = wts[select][subsample]

				if refdata is not None:
					select_data = (refdata.res["snr"]>slow) & (refdata.res["snr"]<shigh) & (refdata.res["mean_rgpp_rp"]>rlow) & (refdata.res["mean_rgpp_rp"]<rhigh) 
					kl = lk.kullback_leibler(eh[select],edat[select_data], show=False)
				else:
					kl=0.0

				# Raise an error if there are too few simulated galaxies in a given bin
				if ngal < 60:
					print "Warning: <100 galaxies in bin %d, %d (ngal=%d)"%(i,j, ngal)
					empty=False
				if ngal==0:
					print "Warning: no galaxies in bin %d, %d "%(i,j)
					empty=True


	
				vrgp_min = rgp_edges[i]
				vsnr_min = snr_edges[j]
				vrgp_max = rgp_edges[i+1]
				vsnr_max = snr_edges[j+1]

				if ngal==0:
					print "Warning: no galaxies in bin %d, %d"%(i,j)
					list_bias.append([j, i, ngal, 0, vrgp_min, vrgp_max, vsnr_min, vsnr_max, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.	, 0., 0.])
					continue

				
				b = di.get_bias(tr[select][subsample], data[select][subsample], weights=bin_wts, apply_calibration=apply_calibration, nbins=10, ellipticity_name=ellipticity_name, binning="equal_number", names=["m","c","m11","m22","c11","c22"], silent=True)
				a = di.get_alpha(data[select][subsample], data[select][subsample], weights=bin_wts, nbins=10, xlim=(-0.03, 0.02), binning="equal_number", names=["alpha", "alpha11", "alpha22"], silent=True, use_weights=False)
				measurement_weight = di.compute_im3shape_weight(data["e1"][select][subsample])

				list_bias.append([j, i, ngal, 10**vrgp_min, 10**vrgp_max, 10**vsnr_min, 10**vsnr_max, measurement_weight, kl, b["m"][0], b["m"][1], b["c"][0], b["c"][1], b["m11"][0], b["m11"][1], b["m22"][0], b["m22"][1], b["c11"][0], b["c11"][1], b["c22"][0], b["c22"][1], a["alpha"][0], a["alpha"][1], a["alpha11"][0], a["alpha11"][1], a["alpha22"][0], a["alpha22"][1] ])

		lab=["j","i","ngal","rgp_lower","rgp_upper","snr_lower","snr_upper","weight","kl","m","err_m","c","err_c","m1","err_m1","m2","err_m2","c1","err_c1","c2","err_c2","alpha","err_alpha","alpha11","err_alpha11","alpha22","err_alpha22"]
		dt = {'names': lab, 'formats': ['i4']*3 + ['f8']*24 }
		arr_bias = np.core.records.fromarrays(np.array(list_bias).transpose(), dtype=dt)

		if table_name is None:
			filename_table_bias = 'bias_table-%s-selection_section%d%s.fits'%(fit, split_half,"_calibrated"*apply_calibration)
		else:
			filename_table_bias = table_name

		import pyfits
		pyfits.writeto(filename_table_bias,arr_bias,clobber=True)
		print 'saved %s'%filename_table_bias


	def diff_table(self,  weights=None, compare_with_subset=None, fit="bord", table_name=None, sbins=10, rbins=5, binning="equal_number", rlim=(1.13,3), slim=(12,200)):
		print 'measuring bias'

		import tools.likelihoods as lk

		data = self.res
		print "using the full catalogue (%d objects)"%(data.size)
		if hasattr(self, "truth"):
			tr = self.truth

		if fit.lower()!="bord":
			print "Using %s only galaxies"%fit
			val = int(fit.lower()=="bulge")
			sel_fit = data["is_bulge"]==val
		else:
			sel_fit = np.ones_like(data).astype(bool)

		data = data[sel_fit]
		tr = tr[sel_fit]
		wts = weights[sel_fit]

		sel_lim = (data["snr"]>slim[0]) & (data["snr"]<slim[1]) & (data["mean_rgpp_rp"]>rlim[0]) & (data["mean_rgpp_rp"]<rlim[1])
		data = data[sel_lim]
		tr = tr[sel_lim]
		wts = wts[sel_lim]
		if compare_with_subset is not None:
			compare_with_subset = compare_with_subset[sel_fit][sel_lim]


		if isinstance(binning,str) : 
			if binning.lower()=="uniform":
				snr_edges = np.logspace(np.log10(slim[0]),np.log10(slim[1]),sbins+1)
				rgp_edges = np.linspace(rlim[0],rlim[1],rbins+1)
			elif binning.lower()=="equal_number":
				snr_edges = di.find_bin_edges(np.log10(data["snr"]), sbins)
				rgp_edges = di.find_bin_edges(np.log10(data["mean_rgpp_rp"]), rbins)

			override_bin_edges = False
			nrbin = rbins
			nsbin = sbins
		elif isinstance(binning, tuple):
			snr_edges = binning[1]
			rgp_edges_low = np.unique(np.array(binning[0]).T[0])
			rgp_edges = np.hstack(( np.unique(np.array(binning[0]).T[0]), np.atleast_1d(np.unique(np.array(binning[0]).T[1])[-1]) ))
			nrbin = np.unique(np.array(binning[0]).T[0]).size
			override_bin_edges = True

		list_bias = []
		bias_grid=[]

		eh =  np.sqrt(data["e1"]*data["e1"] + data["e2"]*data["e2"] )


		if not override_bin_edges:
			print "Will do dynamic binning in SNR"

		for i in xrange(nrbin):
			if not override_bin_edges:	
				snr_samp = data["snr"][(np.log10(data['mean_rgpp_rp']) > rgp_edges[i]) & (np.log10(data['mean_rgpp_rp']) < rgp_edges[i+1])]
				snr_edges=di.find_bin_edges(np.log10(snr_samp), sbins)
			else:
				select_edges = (np.array(binning[0]).T[0]==rgp_edges[i])
				snr_edges_low = np.array(binning[1]).T[0][select_edges] 
				snr_edges_high = np.array(binning[1]).T[1][select_edges]
				snr_edges = np.hstack((snr_edges_low, np.atleast_1d(snr_edges_high[-1])))
				sbins = snr_edges_high.size
			for j in xrange(sbins):
				empty=False
				

				# Select in bins of snr and size
				# This is horrible. I know it, and I'm sorry.
				if not override_bin_edges:
					slow,shigh = 10**snr_edges[j], 10**snr_edges[j+1]
					rlow,rhigh = 10**rgp_edges[i], 10**rgp_edges[i+1]
					select = (data['snr'] > 10**snr_edges[j]) & (data['snr'] < 10**snr_edges[j+1]) & (data['mean_rgpp_rp'] > 10**rgp_edges[i]) & (data['mean_rgpp_rp'] < 10**rgp_edges[i+1])
					print "bin %d %d  snr = [%2.3f-%2.3f] rgpp/rp = [%2.3f-%2.3f]"%(j, i, 10**snr_edges[j], 10**snr_edges[j+1], 10**rgp_edges[i], 10**rgp_edges[i+1] )
				else:
					slow,shigh = 10**snr_edges_low[j], 10**snr_edges_high[j]
					rlow,rhigh = 10**binning[0][i][0], 10**binning[0][i][1]
					select = (data['snr'] > 10**snr_edges_low[j]) & (data['snr'] < 10**snr_edges_high[j]) & (data['mean_rgpp_rp'] > 10**binning[0][i][0]) & (data['mean_rgpp_rp'] < 10**binning[0][i][1])
					print "bin %d %d  snr = [%2.3f-%2.3f] rgpp/rp = [%2.3f-%2.3f]"%(j, i, 10**snr_edges_low[j], 10**snr_edges_high[j], 10**binning[0][i][0], 10**binning[0][i][1] )
				ngal = np.nonzero(select.astype(int))[0].size

				# Raise an error if there are too few simulated galaxies in a given bin
				if ngal < 60:
					print "Warning: <100 galaxies in bin %d, %d (ngal=%d)"%(i,j, ngal)
					empty=False
				if ngal==0:
					print "Warning: no galaxies in bin %d, %d "%(i,j)
					empty=True
	
				vrgp_min = rgp_edges[i]
				vsnr_min = snr_edges[j]
				vrgp_max = rgp_edges[i+1]
				vsnr_max = snr_edges[j+1]

				if ngal==0:
					print "Warning: no galaxies in bin %d, %d"%(i,j)
					list_bias.append([j, i, ngal, 0, vrgp_min, vrgp_max, vsnr_min, vsnr_max, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.	, 0., 0.])
					continue

				
				b = di.get_bias(tr[select], data[select], nbins=10, binning="equal_number", names=["m","c","m11","m22","c11","c22"], silent=True)
				b0 = di.get_bias(tr[select & compare_with_subset], data[select & compare_with_subset], nbins=10, binning="equal_number", names=["m","c","m11","m22","c11","c22"], silent=True)
				a = di.get_alpha(data[select], data[select], nbins=10, xlim=(-0.03, 0.02), binning="equal_number", names=["alpha", "alpha11", "alpha22"], silent=True, use_weights=False)
				a0 = di.get_alpha(data[select & compare_with_subset], data[select & compare_with_subset], nbins=10, xlim=(-0.03, 0.02), binning="equal_number", names=["alpha", "alpha11", "alpha22"], silent=True, use_weights=False)
				measurement_weight = di.compute_im3shape_weight(data["e1"][select])

				kl = lk.kullback_leibler(eh[select],eh[select & compare_with_subset], show=False)

				print "m = [%3.3f %3.3f] %3.4f"%(b["m"][0], b0["m"][0], kl)


				list_bias.append([j, i, ngal, 10**vrgp_min, 10**vrgp_max, 10**vsnr_min, 10**vsnr_max, measurement_weight, kl, b["m"][0], b["m"][1], b0["m"][0], b0["m"][1], a["alpha"][0], a["alpha"][1], a0["alpha"][0], a0["alpha"][1]])

		lab=["j","i","ngal","rgp_lower","rgp_upper","snr_lower","snr_upper","weight","kl","m","err_m","m_subsamp","err_m_subsamp","alpha","err_alpha","alpha_subsamp","err_alpha_subsamp"]
		dt = {'names': lab, 'formats': ['i4']*3 + ['f8']*14 }
		arr_bias = np.core.records.fromarrays(np.array(list_bias).transpose(), dtype=dt)

		if table_name is None:
			filename_table_bias = 'bias_table-%s-selection_section%d%s.fits'%(fit, split_half,"_calibrated"*apply_calibration)
		else:
			filename_table_bias = table_name

		import pyfits
		pyfits.writeto(filename_table_bias,arr_bias,clobber=True)
		print 'saved %s'%filename_table_bias

	def setup_grid_tree(self):
		self.snr = (self.bias_grid["snr_lower"]+self.bias_grid["snr_upper"])/2
		self.rgpp = (self.bias_grid["rgp_lower"]+self.bias_grid["rgp_upper"])/2

		self.fx = np.log10(self.snr).max()
		self.fy = np.log10(self.rgpp).max()

		print "Constructing tree..."
		from scipy.spatial import KDTree as kd
		xy = np.vstack((np.log10(self.snr)/self.fx, np.log10(self.rgpp)/self.fy))
		self.tree = kd(xy.T)

	def get_nearest_node(self):
		if not hasattr(self, "tree"):
			self.setup_grid_tree()

		ngal = self.res["snr"].size

		print "Querying tree..."
		if ngal<1e7:
			xy_target = np.vstack((np.log10(self.res["snr"])/self.fx, np.log10(self.res["mean_rgpp_rp"])/self.fy ))
			self.ids, self.result = self.tree.query(xy_target.T, k=1)
		else:
			print "Splitting array"
			nchunk = ngal/10
			self.ids=[]
			self.result=[]

			for i in xrange(10):
				print i
				xy_target = np.vstack((np.log10(self.res["snr"][i*nchunk:(i+1)*nchunk])/self.fx, np.log10(self.res["mean_rgpp_rp"][i*nchunk:(i+1)*nchunk])/self.fy ))
				ids, result = self.tree.query(xy_target.T, k=1)
				self.ids.append(ids)
				self.result.append(result)
			if ngal!=10*nchunk:
				xy_target = np.vstack((np.log10(self.res["snr"][(i+1)*nchunk:])/self.fx, np.log10(self.res["mean_rgpp_rp"][(i+1)*nchunk:])/self.fy ))
				ids, result = self.tree.query(xy_target.T, k=1)
				self.ids.append(ids)
				self.result.append(result)

	def fit_rbf(self, table="/home/samuroff/bias_table-bord-selection_section1.fits", smoothing=1):
		bt = fi.FITS(table)[1].read()
		self.bt = bt
		self.snr = (bt["snr_lower"]+bt["snr_upper"])/2
		self.rgpp = (bt["rgp_lower"]+bt["rgp_upper"])/2

		print "Fitting to %s"%table

		# Set up the RBF interpolation
		self.fx = np.log10(self.snr).max()
		self.fy = np.log10(self.rgpp).max()
		self.rbf_interp_m = Rbf(np.log10(self.snr)/self.fx, np.log10(self.rgpp)/self.fy,bt["m"], smooth=smoothing, function="multiquadric")
		self.rbf_interp_a = Rbf(np.log10(self.snr)/self.fx, np.log10(self.rgpp)/self.fy,bt["alpha"], smooth=smoothing, function="multiquadric")

	def do_rbf_interpolation(self, bias, cat):

		# If the arrays here  are sufficiently large we may need to apply the interpolation to
		# the catalogue in parts
		interpolators = {"m":self.rbf_interp_m, "a":self.rbf_interp_a}
		if cat.size>10e6:
			print "Splitting array"
			ngal = cat["snr"].size
			nchunk = ngal/5
			interpolated_bias=[]

			for i in xrange(5):
				print i
				interpolated_bias.append( interpolators[bias]( np.log10(cat["snr"][i*nchunk:(i+1)*nchunk])/self.fx, np.log10(cat["mean_rgpp_rp"][i*nchunk:(i+1)*nchunk])/self.fy ) )
			if ngal!=5*nchunk:
				interpolated_bias.append( interpolators[bias]( np.log10(cat["snr"][(i+1)*nchunk:])/self.fx, np.log10(cat["mean_rgpp_rp"][(i+1)*nchunk:])/self.fy ) )

			return np.concatenate(interpolated_bias)

		else:
			return interpolators[bias](np.log10(cat["snr"])/self.fx, np.log10(cat["mean_rgpp_rp"])/self.fy)
		
	def export(self, filename, bias_columns_only=True):
		print "Will write to %s"%filename,
		out = fi.FITS(filename,"rw")

		if bias_columns_only:
			hdu_name = "i3s_calibration_col"
			outcat = np.empty(self.res.size, dtype=[("coadd_objects_id", int), ("m", float), ("c1", float), ("c2", float), ("a", float)])
			for col in outcat.dtype.names:
				outcat[col] = self.res[col]
		else:
			hdu_name = i3s_data
			outcat = self.res

		out.write(outcat)
		out[-1].write_key("EXTNAME", hdu_name)
		out.close()
		print "done"
			
	def fit(self, bias="m", table="/home/samuroff/bias_table-bord-selection_section1.fits", tag=""):
		bt = fi.FITS(table)[1].read()
		snr = (np.unique(bt["snr_lower"])+np.unique(bt["snr_upper"]))/2
		rgpp = (np.unique(bt["rgp_lower"])+np.unique(bt["rgp_upper"]))/2

		mfit = opt.curve_fit(eval_m, np.array([ (bt["snr_lower"]+bt["snr_upper"])/2,(bt["rgp_lower"]+bt["rgp_upper"])/2 ]).T, bt["m"],p0=wm_tk) 
		afit = opt.curve_fit(eval_a, np.array([ (bt["snr_lower"]+bt["snr_upper"])/2,(bt["rgp_lower"]+bt["rgp_upper"])/2 ]).T, bt["alpha"],p0=wm_tk) 
		if bias=="m":
			par = mfit[0]
		elif bias=="a":
			par = afit[0]

		setattr(self, "optimised_coefficients_%s"%bias, par)
		x = np.array([ (bt["snr_lower"]+bt["snr_upper"])/2,(bt["rgp_lower"]+bt["rgp_upper"])/2 ]).T

		com2="a0 "+", a%d "*par[1:].size +"=tuple(par)"
		com2=com2%tuple(np.linspace(1,17,17))
		exec com2

		com="eval_%s(x"%bias + ", a%d "*par.size+")"
		com=com%tuple(np.linspace(0,17,18))
		exec "bias=%s"%com

	def merge_bias_cols(self, nbc_disc, nbc_bulge, split_half=2, bias="m", mask=None):

		if split_half==0:
			res = getattr(self, "res")
			resb = getattr(nbc_bulge, "res")
			resd = getattr(nbc_disc, "res")
		else:
			res = getattr(self, "res%d"%split_half)
			resb = getattr(nbc_bulge, "res%d"%split_half)
			resd = getattr(nbc_disc, "res%d"%split_half)

		if mask is None:
			mask = np.ones_like(res['e1']).astype(bool)

		print "Checking columns"

		if not ("m" in res.dtype.names):
			res = arr.add_col(res, "m", np.empty(res['e1'].size), verbose=False, dtype=np.float64)
			res = arr.add_col(res, "c1", np.empty(res['e1'].size), verbose=False, dtype=np.float64)
			res = arr.add_col(res, "c2", np.empty(res['e1'].size), verbose=False, dtype=np.float64)

		bulge0 = res["is_bulge"].astype(bool)
		bulge1 = resb["is_bulge"].astype(bool)
		print "column : %s, bulge : %d/%d, disc : %d/%d"%(bias, res[mask & bulge0].size, res[mask].size, res[mask & np.invert(bulge0)].size, res[mask].size)
		
		res[bias][mask & bulge0] = resb[bias][bulge1]
		
		res[bias][mask & np.invert(bulge0)] = resd[bias][np.invert(bulge1)]

		return res

	def combine_bd(self, nbc_disc, nbc_bulge, split_half=2, names=["m", "c1", "c2"], mask=None):
		print "Will combine bulge and disc calibration fits."

		suffix = "" + "1"*(split_half==1) + "2"*(split_half==2)

		for bias_name in names:
			setattr(self, "res%s"%suffix, self.merge_bias_cols(nbc_disc,  nbc_bulge, split_half=split_half, bias=bias_name, mask=mask) )
			gc.collect()

	def get_combined_calibration(self, nbc_disc, nbc_bulge, split_half=2, names=["m", "c1", "c2"], mask=None):
		print "Will combine bulge and disc calibration fits."
		if split_half==0:
			if mask is None: mask=np.ones_like(self.res['e1']).astype(bool)
			for bias in names:
				if "m" not in self.res.dtype.names:
					self.res = arr.add_col(self.res, "m", np.zeros_like(self.res['e1']))
					self.res = arr.add_col(self.res, "c1", np.zeros_like(self.res['e1']))
					self.res = arr.add_col(self.res, "c2", np.zeros_like(self.res['e1']))
				bulge = self.res[mask]["is_bulge"].astype(bool)
				print "column : %s, bulge : %d/%d, disc : %d/%d"%(bias, self.res[mask][bulge].size, self.res[mask].size, self.res[mask][np.invert(bulge)].size, self.res[mask].size)
				try:
					self.res[mask][bias][bulge] = nbc_bulge.res[bias][bulge]
				except:
					import pdb ; pdb.set_trace()
				self.res[mask][bias][np.invert(bulge)] = nbc_disc.res[bias][np.invert(bulge)]

		else:
			
			com ="""
for i, bias in enumerate(names):
	if mask is None: mask=np.ones_like(self.res['e1']).astype(bool)
	bulge = self.res[mask]['is_bulge'].astype(bool)
	if i==0: print 'bulge :', self.res[mask][bulge].size, 'disc : ', self.res[mask][np.invert(bulge)].size, 'total : ', self.res[mask].size
	self.res = arr.add_col(self.res, bias, np.zeros_like(self.res['e1']))

	print 'column : ', bias
				
	self.res[mask][bias][bulge] = nbc_bulge.res[bias][bulge]
	self.res[mask][bias][np.invert(bulge)] = nbc_disc.res[bias][np.invert(bulge)]""".replace("res", "res%d"%split_half)
			exec(com)
		print "done"



def show_table(table_name,ls="none", fmt="o", legend=False, name="m", do_half=0):
	bt = fi.FITS(table_name)[1].read()
	rgpp = (np.unique(bt["rgp_lower"])+np.unique(bt["rgp_upper"]))/2
	nbins = rgpp.size

	plt.xscale("log")
	colours=["purple", "forestgreen", "steelblue", "pink", "darkred", "midnightblue", "gray", "sienna", "olive", "darkviolet","cyan", "red", "k", "deepskyblue", "yellow", "darkseagreen", "darksalmon"]
	pts = ["o", "D", "x", "^", ">", "<", "1", "s", "*", "+", ".", "o", "D", "x", "^", ">", "<", "1", "s", "*", "+", "."]
	for i,r in enumerate(rgpp):
		sel = (bt["i"]==i)
		snr = 10** ((np.log10(bt["snr_lower"][sel]) + np.log10(bt["snr_upper"][sel]))/2)

		if do_half==1 and i>nbins/2:
			continue
		elif do_half==2 and i<nbins/2:
			continue
		if legend:
			try:
				plt.errorbar(snr, bt["%s"%name][i*snr.size:(i*snr.size)+snr.size], bt["err_%s"%name][i*snr.size:(i*snr.size)+snr.size], color=colours[i], ls=ls, fmt=pts[i], lw=2.5, label="$R_{gpp}/R_p = %1.2f-%1.2f$"%(np.unique(bt["rgp_lower"])[i],np.unique(bt["rgp_upper"])[i]))
			except:
				import pdb ; pdb.set_trace()
		else:
			plt.errorbar(snr, bt["%s"%name][i*snr.size:(i*snr.size)+snr.size], bt["err_%s"%name][i*snr.size:(i*snr.size)+snr.size], color=colours[i], ls=ls, fmt=pts[i], lw=2.5)

	plt.xlim(10,250)
	plt.axhline(0, lw=2, color="k")
	
	plt.xlabel("Signal-to-Noise $SNR_w$")
	if name=="m":
		plt.ylim(-0.85,0.05)
		plt.ylabel("Multiplicative Bias $m \equiv (m_1 + m_2)/2$")
	elif name=="alpha":
		plt.ylabel(r"PSF Leakage $\alpha \equiv (\alpha _1 + \alpha _2)/2$")
		plt.ylim(-0.5,2)



	plt.legend(loc="lower right")


wm_tk = np.array([-1.05e-03,1.47e-06,8.10e-05,6.73e-06,0.00e+00,0.00e+00,0.00e+00,6.56e-05,-6.16e-06,-2.09e-05,-7.63e-03,-1.37e-03,-1.08e-04,1.63e-03,5.37e-03,1.63e-04,2.28e-03,-2.73e-11])
wa_tk = np.array([-1.67283612e-04, 1.09715332e-06, 5.95801408e-05, 6.39015150e-07, 2.97121531e-08, -3.60228146e-10, 4.73608639e-09, 4.05443791e-05, -3.52379986e-06, -1.95627195e-05, 8.97549111e-04, -3.23420375e-04, -1.91942923e-06, -6.57971727e-12, -1.41012000e-09, 1.61504257e-15, 2.36381064e-11, -1.76498862e-12])

def basis_m(snr_rgpp):
	snr1=snr_rgpp[:,0]/100.
	rgp1=(snr_rgpp[:,1]-1.)/10.

	func=np.zeros((len(snr1),18))
	func[:,0]=1/snr1**2*1/rgp1**2
	func[:,1]=1/snr1**3*1/rgp1**3
	func[:,2]=1/snr1**3*1/rgp1**2
	func[:,3]=1/snr1**2*1/rgp1**3
	func[:,4]=1/snr1**4*1/rgp1**3
	func[:,5]=1/snr1**4*1/rgp1**4
	func[:,6]=1/snr1**3*1/rgp1**4
	func[:,7]=1/snr1**2.5*1/rgp1**2.5
	func[:,8]=1/snr1**2.5*1/rgp1**3
	func[:,9]=1/snr1**3*1/rgp1**2.5
	func[:,10]=1/snr1**1.5*1/rgp1**1.5
	func[:,11]=1/snr1**2.5*1/rgp1**1.5
	func[:,12]=1/snr1**1.5*1/rgp1**2.5
	func[:,13]=1/snr1**1.5*1/rgp1**2
	func[:,14]=1/snr1**2*1/rgp1**1.5
	func[:,15]=1/snr1**1.25*1/rgp1**1.75
	func[:,16]=1/snr1**1.75*1/rgp1**1.25
	func[:,17]=1/snr1**4*1/rgp1**4

	return func

def eval_m(x, a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17):
	basis = basis_m(x)
	return np.dot(basis, np.array([a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17]) )

def basis_a(snr_rgpp):
	snr1=snr_rgpp[:,0]/100.
	rgp1=(snr_rgpp[:,1]-1.)/10.

	func=np.zeros((len(snr1),18))
	func[:,0] =1/snr1**2*1/rgp1**2
	func[:,1] =1/snr1**3*1/rgp1**3
	func[:,2] =1/snr1**3*1/rgp1**2
	func[:,3] =1/snr1**2*1/rgp1**3
	func[:,4] =1/snr1**4*1/rgp1**3
	func[:,5] =1/snr1**4*1/rgp1**4
	func[:,6] =1/snr1**3*1/rgp1**4
	func[:,7] =1/snr1**2.5*1/rgp1**2.5
	func[:,8] =1/snr1**2.5*1/rgp1**3
	func[:,9] =1/snr1**3* 1/rgp1**2.5
	func[:,10]=1/snr1**1.5*1/rgp1**1.5
	func[:,11]=1/snr1**2.5*1/rgp1**1.5
	func[:,12]=1/snr1**1.5*1/rgp1**2.5
	func[:,13]=1/snr1**3*1/rgp1**5
	func[:,14]=1/snr1**5*1/rgp1**3
	func[:,15]=1/snr1**5*1/rgp1**5
	func[:,16]=1/snr1**5*1/rgp1**4
	func[:,17]=1/snr1**4*1/rgp1**5

	return func

def eval_a(x, a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17):
	basis = basis_a(x)
	return np.dot(basis, np.array([a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17]) )

def fit_model(xt, yt, zt, st, ng):
	from tktools import fitting

	xx=np.concatenate([xt[:,None],yt[:,None]],axis=1)
	w,w_cov = get_optimal_w(xx, zt, st, n_gals=ng, expand_basis=basis_m)
	b,sigma = fitting.predict(xx,w,w_cov,expand=basis_m)

	r = b-zt

	print 'expected error=%2.4f' % (np.sum(r*ng)/float(np.sum(ng)))

	return w, w_cov

def other_half(split):
	if split==1:
		return 2
	elif split==2:
		return 1
	else: 
		return None

def get_optimal_w(snr_mid, m, s, expand_basis, n_gals=1):
	list_chi2 = []
	list_res_sig = []
	list_w = []
	list_w_cov = []
	eps_grid=np.linspace(-15,15,1000)

	for eps in eps_grid:
		w, w_cov = tktools.fitting.fit(snr_mid, m, s,expand=expand_basis,eps=10**eps)
		p_mid, sigma = tktools.fitting.predict(snr_mid,w,w_cov,expand=expand_basis)
		chi2 =  np.mean( (((m-p_mid)**2)/(s**2)) )  - np.sum( n_gals*(m-p_mid) )/ float(np.sum(n_gals))
		# chi2 =  np.mean( (((m-p_mid)**2)/(s**2)) ) 
		res_sig = np.sum( n_gals*(m-p_mid) )/ float(np.sum(n_gals))
		# print 'optimal w: eps=%2.2e chi2=%2.2e mean_res/sig=%2.4e' % (10**(eps),chi2,res_sig)
		list_chi2.append(chi2)
		list_w.append(w)
		list_w_cov.append(w_cov)
		list_res_sig.append(res_sig)

		arr_chi2 = np.array(list_chi2)

		select = np.argmin(np.abs(arr_chi2-1))

		w = list_w[select]
		w_cov = list_w_cov[select]
		eps = eps_grid[select]
		chi2 = list_chi2[select]
		res_sig = list_res_sig[select]

		print 'final optimal w: eps=%2.2e chi2=%2.2f res_sig=%2.5f' % (10**(eps),chi2,res_sig)

		return w, w_cov
