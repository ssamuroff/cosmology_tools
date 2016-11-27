import numpy as np
import scipy as sp
import astropy.io.fits as pyfits
import os, pdb
import tools.diagnostics as di
import tools.arrays as arr
import tools.plots as plots
import tools.shapes as sh
import tktools
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

	def get_split_data(self, cat):
		self.res = cat.res
		sel = np.random.randint(0,2,cat.res.size).astype(bool)
		self.res1 = self.res[sel]
		self.res2 = self.res[np.invert(sel)]

		if hasattr(cat,"truth"):
			self.truth = cat.truth
			self.truth1 = cat.truth[sel]
			self.truth2 = cat.truth[np.invert(sel)]

	def apply(self, names=["m","a"], split_half=2, use_rbf=False):

		if use_rbf:
			print "Using RBF interpolation"
		else:
			print "Using polynomial fit"

		# Choose input data
		if split_half>0:
			exec "catalogue_to_calibrate = self.res%d"%split_half
		else:
			catalogue_to_calibrate = self.res

		# Apply the best-fit calibration coefficients

		for bias_name in names:
			print "creating column %s"%bias_name
			if not use_rbf:
				com2 = "a0 "+", a%d "*self.optimised_coefficients_m[1:].size +"=tuple(self.optimised_coefficients_%s)"%bias_name
				com2 = com2%tuple(np.linspace(1,17,17))
				exec com2

				com="eval_%s(np.array([catalogue_to_calibrate['snr'],catalogue_to_calibrate['mean_rgpp_rp']]).T"%bias_name + ", a%d "*self.optimised_coefficients_m.size+")"
				com=com%tuple(np.linspace(0,17,18))
				exec "bias=%s"%com
			else:
				try:
					bias = self.do_rbf_interpolation(bias_name, catalogue_to_calibrate)
				except:
					import pdb ; pdb.set_trace()

			if split_half>0:
				exec "self.res%d = arr.add_col(self.res2, '%s', bias)"%(split_half, bias_name)
				if bias_name=="a":
					exec "self.res%d = arr.add_col(self.res2, 'c1', bias*self.res%d['mean_psf_e1_sky'])"%(split_half, split_half)
					exec "self.res%d = arr.add_col(self.res2, 'c2', bias*self.res%d['mean_psf_e2_sky'])"%(split_half, split_half)
			else:
				self.res = arr.add_col(self.res, bias_name, bias)
				if bias_name=="a":
					self.res = arr.add_col(self.res, "c1", bias*self.res["mean_psf_e1_sky"])
					self.res = arr.add_col(self.res, "c2", bias*self.res["mean_psf_e2_sky"])

			print "Finished calibration"

	def compute(self, split_half=0, fit="bord", apply_calibration=False, table_name=None, ellipticity_name="e", sbins=10, rbins=5, binning="equal_number",rlim=(1,3), slim=(10,1000)):
		print 'measuring bias'

		if split_half>0:
			exec "data = self.res%d"%split_half
			print "using half %d of the catalogue (%d objects)"%(split_half,data.size)
			if hasattr(self, "truth"):
				exec "tr = self.truth%d"%split_half
		else:
			data = self.res[sel_fit]
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

		sel_lim = (data["snr"]>slim[0]) & (data["snr"]<slim[1]) & (data["mean_rgpp_rp"]>rlim[0]) & (data["mean_rgpp_rp"]<rlim[1])
		data = data[sel_lim]
		tr = tr[sel_lim]

		if isinstance(binning,str) : 
			if binning.lower()=="uniform":
				snr_edges = np.logspace(np.log10(slim[0]),np.log10(slim[1]),sbins+1)
				rgp_edges = np.linspace(rlim[0],rlim[1],rbins+1)
			elif binning.lower()=="equal_number":
				snr_edges = di.find_bin_edges(np.log10(data["snr"]), sbins)
				rgp_edges = di.find_bin_edges(np.log10(data["mean_rgpp_rp"]), rbins)

		snr_centres = (snr_edges[1:]+snr_edges[:-1])/2.0
		rgp_centres = (rgp_edges[1:]+rgp_edges[:-1])/2.0

		list_bias = []
		bias_grid=[]

		b = di.get_bias(tr, data, nbins=5, apply_calibration=apply_calibration, ellipticity_name=ellipticity_name, binning="equal_number", names=["m","c","m11","m22","c11","c22"], silent=True)
		print "Global biases:"
		print "m11 : ", b["m11"]
		print "m22 : ", b["m22"]
		print "m : ", b["m"]
		print "c11 : ", b["c11"]
		print "c22 : ", b["c22"]
		print "c : ", b["c"]

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
					list_bias.append([j, i, ngal, 0, vrgp_min, vrgp_max, vsnr_min, vsnr_max, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.	, 0., 0.])
					continue

				

				filename_str = 'snr%2.2f.rgpp%2.2f' % (10**vsnr_mid,10**vrgp_mid)
				b = di.get_bias(tr[select], data[select], apply_calibration=apply_calibration, nbins=5, ellipticity_name=ellipticity_name, binning="equal_number", names=["m","c","m11","m22","c11","c22"], silent=True)
				a = di.get_alpha(data[select], data[select], nbins=5, xlim=(-0.015, 0.02), binning="equal_number", names=["alpha", "alpha11", "alpha22"], silent=True, use_weights=False)

				list_bias.append([j, i, ngal, 10**vrgp_min, 10**vrgp_max, 10**vsnr_min, 10**vsnr_max, b["m"][0], b["m"][1], b["c"][0], b["c"][1], b["m11"][0], b["m11"][1], b["m22"][0], b["m22"][1], b["c11"][0], b["c11"][1], b["c22"][0], b["c22"][1], a["alpha"][0], a["alpha"][1], a["alpha11"][0], a["alpha11"][1], a["alpha22"][0], a["alpha22"][1] ])

		lab=["j","i","ngal","rgp_lower","rgp_upper","snr_lower","snr_upper","m","err_m","c","err_c","m1","err_m1","m2","err_m2","c1","err_c1","c2",	"err_c2","alpha","err_alpha","alpha11","err_alpha11","alpha22","err_alpha22"]
		dt = {'names': lab, 'formats': ['i4']*3 + ['f8']*22 }
		arr_bias = np.core.records.fromarrays(np.array(list_bias).transpose(), dtype=dt)

		if table_name is None:
			filename_table_bias = 'bias_table-%s-selection_section%d%s.fits'%(fit, split_half,"_calibrated"*apply_calibration)
		else:
			filename_table_bias = table_name

		import pyfits
		pyfits.writeto(filename_table_bias,arr_bias,clobber=True)
		print 'saved %s'%filename_table_bias

	def fit_rbf(self, table="/home/samuroff/bias_table-bord-selection_section1.fits"):
		bt = fi.FITS(table)[1].read()
		snr = (bt["snr_lower"]+bt["snr_upper"])/2
		rgpp = (bt["rgp_lower"]+bt["rgp_upper"])/2

		# Set up the RBF interpolation
		self.fx = np.log10(snr).max()
		self.fy = np.log10(rgpp).max()
		self.rbf_interp_m = Rbf(np.log10(snr)/self.fx, np.log10(rgpp)/self.fy,bt["m"], smooth=3, function="multiquadric")
		self.rbf_interp_a = Rbf(np.log10(snr)/self.fx, np.log10(rgpp)/self.fy,bt["alpha"], smooth=3, function="multiquadric")

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

	def export(self, filename):
		print "Will write to %s"%filename,
		out = fi.FITS(filename,"rw")
		out.write(self.res)
		out.write_header("EXTNAME", "i3s_data")
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


	def get_combined_calibration(self, nbc_disc, nbc_bulge, split_half=2, names=["m", "c1", "c2"]):
		print "Will combine bulge and disc calibration fits."
		if split_half==0:
			for bias in names:
				self.res = arr.add_col(self.res, bias, np.zeros_like(self.res['e1']))
				bulge = self.res["is_bulge"].astype(bool)
				print "column : %s, bulge : %d/%d, disc : %d/%d"%(bias, self.res[bulge].size, self.res.size, self.res[np.invert(bulge)].size, self.res.size)
				try:
					self.res[bias][bulge] = nbc_bulge.res[bias][bulge]
				except:
					import pdb ; pdb.set_trace()
				self.res[bias][np.invert(bulge)] = nbc_disc.res[bias][np.invert(bulge)]

		else:
			
			com ="""
for i, bias in enumerate(names):
	bulge = self.res['is_bulge'].astype(bool)
	if i==0: print 'bulge :', self.res[bulge].size, 'disc : ', self.res[bulge].size, 'total : ', self.res.size
	self.res = arr.add_col(self.res, bias, np.zeros_like(self.res['e1']))

	print 'column : ', bias
				
	self.res[bias][bulge] = nbc_bulge.res[bias][bulge]
	self.res[bias][np.invert(bulge)] = nbc_disc.res[bias][np.invert(bulge)]""".replace("res", "res%d"%split_half)
			exec(com)
		print "done"

def show_table(table_name,ls="none", fmt="o", legend=False, name="m", do_half=0):
	bt = fi.FITS(table_name)[1].read()
	rgpp = (np.unique(bt["rgp_lower"])+np.unique(bt["rgp_upper"]))/2
	nbins = rgpp.size

	plt.xscale("log")
	colours=["purple", "forestgreen", "steelblue", "pink", "darkred", "midnightblue", "gray", "sienna", "olive", "darkviolet"]
	pts = ["o", "D", "x", "^", ">", "<", "1", "s", "*", "+", "."]
	for i,r in enumerate(rgpp):
		sel = (bt["i"]==i)
		snr = 10** ((np.log10(bt["snr_lower"][sel]) + np.log10(bt["snr_upper"][sel]))/2)

		if do_half==1 and i>nbins/2:
			continue
		elif do_half==2 and i<nbins/2:
			continue
		if legend:
			plt.errorbar(snr, bt["%s"%name][i*snr.size:(i*snr.size)+snr.size], bt["err_%s"%name][i*snr.size:(i*snr.size)+snr.size], color=colours[i], ls=ls, fmt=pts[i], lw=2.5, label="$R_{gpp}/R_p = %1.2f-%1.2f$"%(np.unique(bt["rgp_lower"])[i],np.unique(bt["rgp_upper"])[i]))
		else:
			plt.errorbar(snr, bt["%s"%name][i*snr.size:(i*snr.size)+snr.size], bt["err_%s"%name][i*snr.size:(i*snr.size)+snr.size], color=colours[i], ls=ls, fmt=pts[i], lw=2.5)

	plt.xlim(10,300)
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
