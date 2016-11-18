import numpy as np
import scipy as sp
import pyfits
import os, pdb
from samplers import sampler
import tools.emcee as mc
import tools.diagnostics as di
import tools.arrays as arr
import tools.plots as plots
import tktools

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

class nbc(plots.im3shape_results_plots):
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

	def apply(self, names=["m","a"], split_half=2):

		# Choose input data
		if split_half>0:
			exec "catalogue_to_calibrate = self.res%d"%split_half
		else:
			catalogue_to_calibrate = self.res

		# Apply the best-fit calibration coefficients
		x = np.array([catalogue_to_calibrate["snr"],catalogue_to_calibrate["mean_rgpp_rp"]]).T

		for bias_name in names:
			print "creating column %s"%bias_name
			com2 = "a0 "+", a%d "*self.optimised_coefficients_m[1:].size +"=tuple(self.optimised_coefficients_%s)"%bias_name
			com2 = com2%tuple(np.linspace(1,17,17))
			exec com2

			com="eval_%s(x"%bias_name + ", a%d "*self.optimised_coefficients_m.size+")"
			com=com%tuple(np.linspace(0,17,18))
			exec "bias=%s"%com

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



	def compute(self, split_half=0, apply_calibration=False, ellipticity_name="e", sbins=10, rbins=5, binning="equal_number", error_type="bootstrap",rlim=(1,3), slim=(10,1000)):
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

		if hasattr(self, "truth"):
			data = arr.add_col(data,"true_g1",tr["true_g1"])
			data = arr.add_col(data,"true_g2",tr["true_g2"])
			data = arr.add_col(data,"intrinsic_sheared_e1",tr["intrinsic_e1"]+tr["true_g1"])
			data = arr.add_col(data,"intrinsic_sheared_e2",tr["intrinsic_e2"]+tr["true_g2"])

		if isinstance(binning,str) : 
			if binning.lower()=="uniform":
				snr_edges = np.logspace(np.log10(slim[0]),np.log10(slim[1]),sbins+1)
				rgp_edges = np.linspace(rlim[0],rlim[1],rbins+1)
			elif binning.lower()=="equal_number":
				snr_edges=di.find_bin_edges(self.res["snr"][(self.res["snr"]>slim[0]) & (self.res["snr"]<slim[1])], sbins)
				rgp_edges=di.find_bin_edges(self.res["mean_rgpp_rp"][(self.res["mean_rgpp_rp"]>rlim[0]) & (self.res["mean_rgpp_rp"]<rlim[1])], rbins)

		snr_centres = (snr_edges[1:]+snr_edges[:-1])/2.0
		rgp_centres = (rgp_edges[1:]+rgp_edges[:-1])/2.0

		list_bias = []
		bias_grid=[]

		b = di.get_bias(data, nbins=5, apply_calibration=apply_calibration, ellipticity_name=ellipticity_name, binning="equal_number", names=["m","c","m11","m22","c11","c22"], silent=True)
		print "Global biases:"
		print "m11 : ", b["m11"]
		print "m22 : ", b["m22"]
		print "m : ", b["m"]
		print "c11 : ", b["c11"]
		print "c22 : ", b["c22"]
		print "c : ", b["c"]

		for i in xrange(len(rgp_edges)-1):
			for j in xrange(len(snr_edges)-1):
				empty=False
				print "bin %d %d  snr = [%2.3f-%2.3f] rgpp/rp = [%2.3f-%2.3f]"%(j, i, snr_edges[j], snr_edges[j+1], rgp_edges[i], rgp_edges[i+1] )

				# Select in bins of snr and size
				select = (data['snr'] > snr_edges[j]) & (data['snr'] < snr_edges[j+1]) & (data['mean_rgpp_rp'] > rgp_edges[i]) & (data['mean_rgpp_rp'] < rgp_edges[i+1])
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

				filename_str = 'snr%2.2f.rgpp%2.2f' % (vsnr_mid,vrgp_mid)
				b = di.get_bias(data[select], apply_calibration=apply_calibration, nbins=5, ellipticity_name=ellipticity_name, binning="equal_number", names=["m","c","m11","m22","c11","c22"], silent=True)
				a = di.get_alpha(data[select], nbins=5, binning="equal_number", names=["alpha", "alpha11", "alpha22"], silent=True)

				list_bias.append([j, i, ngal, vrgp_min, vrgp_max, vsnr_min, vsnr_max, b["m"][0], b["m"][1], b["c"][0], b["c"][1], b["m11"][0], b["m11"][1], b["m22"][0], b["m22"][1], b["c11"][0], b["c11"][1], b["c22"][0], b["c22"][1], a["alpha"][0], a["alpha"][1], a[	"alpha11"][0], a["alpha11"][1], a["alpha22"][0], a["alpha22"][1] ])

		lab=["j","i","ngal","rgp_lower","rgp_upper","snr_lower","snr_upper","m","err_m","c","err_c","m1","err_m1","m2","err_m2","c1","err_c1","c2",	"err_c2","alpha","err_alpha","alpha11","err_alpha11","alpha22","err_alpha22"]
		dt = {'names': lab, 'formats': ['i4']*3 + ['f8']*22 }
		arr_bias = np.core.records.fromarrays(np.array(list_bias).transpose(), dtype=dt)

		filename_table_bias = 'bias_table-bord-selection_section%d%s.fits'%(split_half,"_calibrated"*apply_calibration)

		import pyfits
		pyfits.writeto(filename_table_bias,arr_bias,clobber=True)
		print 'saved %s'%filename_table_bias

	def fit(self, bias="m", table="/home/samuroff/bias_table-bord-selection_section1.fits", tag="", visual=True):
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

		if visual:
			plt.xscale("log")
			colours=["purple", "forestgreen", "steelblue", "pink", "darkred", "midnightblue"]
			for i,r in enumerate(rgpp):
				plt.plot(x.T[0][i*snr.size:(i*snr.size)+snr.size], bias[i*snr.size:(i*snr.size)+snr.size], color=colours[i], lw=2.5, label="$R_{gpp}/R_p = %1.2f-%1.2f$"%(np.unique(bt["rgp_lower"])[i],np.unique(bt["rgp_upper"])[i]))

			plt.xlim(10,200)
			plt.axhline(0, lw=2)
			plt.ylim(-0.5,0.05)
			plt.legend(loc="lower right")
			plt.savefig("/home/samuroff/shear_pipeline/bias_calibration/plots/v2.2/fits/biasfits_v2.2%s.png"%tag)
			plt.close()

def show_table(table_name,ls="-", legend=False):
	bt = fi.FITS(table_name)[1].read()
	snr = (np.unique(bt["snr_lower"])+np.unique(bt["snr_upper"]))/2
	rgpp = (np.unique(bt["rgp_lower"])+np.unique(bt["rgp_upper"]))/2
	x = np.array([ (bt["snr_lower"]+bt["snr_upper"])/2,(bt["rgp_lower"]+bt["rgp_upper"])/2 ]).T

	plt.xscale("log")
	colours=["purple", "forestgreen", "steelblue", "pink", "darkred", "midnightblue"]
	for i,r in enumerate(rgpp):
		if legend:
			plt.errorbar(x.T[0][i*snr.size:(i*snr.size)+snr.size], bt["m"][i*snr.size:(i*snr.size)+snr.size], bt["err_m"][i*snr.size:(i*snr.size)+snr.size], color=colours[i], ls=ls, lw=2.5, label="$R_{gpp}/R_p = %1.2f-%1.2f$"%(np.unique(bt["rgp_lower"])[i],np.unique(bt["rgp_upper"])[i]))
		else:
			plt.errorbar(x.T[0][i*snr.size:(i*snr.size)+snr.size], bt["m"][i*snr.size:(i*snr.size)+snr.size], bt["err_m"][i*snr.size:(i*snr.size)+snr.size], color=colours[i], ls=ls, lw=2.5)

	plt.xlim(10,300)
	plt.axhline(0, lw=2)
	plt.ylim(-0.65,0.1)
	plt.xlabel("Signal-to-Noise $SNR_w$")
	plt.ylabel("Multiplicative Bias $m \equiv (m_1 + m_2)/2$")
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