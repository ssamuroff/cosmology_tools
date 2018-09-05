import numpy as np
import tools.diagnostics as di
import tools.plots as pl
import fitsio as fi
import os, yaml, argparse, glob
import tools.shapes as s
from tools.im3shape import calibrate_all_tomographic as ct
from tools.im3shape import calibrate_all as ca

import matplotlib.colors
import matplotlib
import pylab as plt
from matplotlib.patches import Ellipse
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['font.size']=16
matplotlib.rcParams['legend.fontsize']=15
matplotlib.rcParams['xtick.major.size'] = 10.0
matplotlib.rcParams['ytick.major.size'] = 10.0

class paper:
	def __init__(self):
		print "tool for remaking paper plots"
		self.config = yaml.load(open("/home/samuroff/local/python/lib/python2.7/site-packages/tools/im3shape/calibration/config/fiducial-y1-unblinded.yaml"))

	def load_data(self, data=False, waxwing=False):
		if not hasattr(self, "hoopoe"):
			self.hoopoe, self.weights, self.y1v2 = ct.setup(True, data, self.config)
		if not hasattr(self, "waxwing"):
			self.waxwing = s.shapecat(res="/share/des/disc8/cambridge/bord-fits/combined_cats/waxwing-500-ohio_A4_A6-infocuts-hv2-i3sv2-catv02.fits", truth="/share/des/disc8/cambridge/bord-fits/combined_cats/waxwing-500-ohio_A4_A6-infocuts-hv2-i3sv2-catv02.fits")

			self.waxwing.load(truth=True)

	def fig1(self, cid=3007825268, tile="DES2111+0043", iobj=None):
		mh=s.meds_wrapper(glob.glob("/share/des/disc8/cambridge/meds/%s-*.fits.fz"%tile)[0], model="disc")
		m0=s.meds_wrapper(glob.glob("/home/samuroff/%s-*.fits.fz"%tile)[0], model="disc")

		if iobj is None:
			iobj=np.argwhere(mh._fits["object_data"].read()["id"]==cid)[0,0]
		else:
			cid = mh._fits["object_data"].read()["id"][iobj]

		res,image,model,wts, transform, vals = mh.i3s(0, coadd_objects_id=cid, return_vals=True, show=False)
		res0,image0,model0,wts0, transform0, vals0 = m0.i3s(0, coadd_objects_id=cid, return_vals=True, show=False)
		
		inputs = mh.get_im3shape_inputs(iobj)
		transforms = inputs.all('transform')
		images = mh.get_cutout_list(iobj)[1:]
		images0 = m0.get_cutout_list(iobj)[1:]
		nexp = len(images)

		offsets=[]
		ellipses=[]

		fig=plt.figure()

		e = res.e1*res.e1 + res.e2*res.e2
		e = np.sqrt(e)
		e0 = res.e1 + 1j*res.e2
		phi = np.angle(e0, deg=True)/2
		q=(1-e)/(1+e)
		b=7
		a=b*q

		A=np.array([np.array([0.,0.]),np.array([0.,0.])])

		iplot=1


		for i, transform in enumerate(transforms):
			print iplot
			A[0,0]=transform.A[0][0]
			A[0,1]=transform.A[0][1]
			A[1,1]=transform.A[1][1]
			A[1,0]=transform.A[1][0]
			xf,yf = np.dot([res.ra_as,res.dec_as], A)

			ellip = Ellipse((int(transform.x0+xf), int(transform.y0+yf)), 2*a, 2*b, phi, edgecolor='purple',facecolor='none')

			ax = fig.add_subplot(3, nexp, iplot)
			ax.imshow(images[i], interpolation="none")
			plt.setp(ax.get_yticklabels(), visible=False)
			plt.setp(ax.get_xticklabels(), visible=False)
			ax.add_artist(ellip)
			ellip.set_facecolor("none")
			ellip.set_edgecolor("purple")
			ellip.set_linewidth(2.5)
			iplot+=1

		e = res0.e1*res0.e1 + res0.e2*res0.e2
		e = np.sqrt(e)
		e0 = res0.e1 + 1j*res0.e2
		phi = np.angle(e0, deg=True)/2
		q=(1-e)/(1+e)
		b=7
		a=b*q

		A=np.array([np.array([0.,0.]),np.array([0.,0.])])


		for i, transform in enumerate(transforms):
			print iplot
			A[0,0]=transform.A[0][0]
			A[0,1]=transform.A[0][1]
			A[1,1]=transform.A[1][1]
			A[1,0]=transform.A[1][0]
			xf,yf = np.dot([res0.ra_as,res0.dec_as], A)

		

			ellip = Ellipse((int(transform.x0+xf), int(transform.y0+yf)), 2*a, 2*b, phi, edgecolor='purple',facecolor='none')

			ax = fig.add_subplot(3, nexp, iplot)
			ax.imshow(images0[i], interpolation="none")
			plt.setp(ax.get_yticklabels(), visible=False)
			plt.setp(ax.get_xticklabels(), visible=False)
			ax.add_artist(ellip)
			ellip.set_facecolor("none")
			ellip.set_edgecolor("purple")
			ellip.set_linewidth(2.5)
			iplot+=1

		plt.subplots_adjust(wspace=0., hspace=0.0)

		plt.savefig("%s/%s"%(".","simple_hoopoe_waxwing_example_%s_obj%d.png"%(tile,cid)))
		plt.close()
		return vals, vals0


	def fig6(self, generate_data_snr=True,generate_data_r=True, alpha=False, generate_plot=True):
	
		if generate_data_snr:
			# blinded
			pth="/share/des/disc8/cambridge/bord-fits/combined_cats/matched_hoopoe-waxwing_nocuts.fits"
			# unblinded
			pth="/share/des/disc8/y1-unblinded-combined_cats/hoopoe-waxwing-overlap-unblinded.fits"
			print pth
			matched_hoopoe =s.shapecat(res=pth)
			matched_waxwing =s.shapecat(res=pth)
			hoopoe =s.shapecat(res=pth)
			waxwing =s.shapecat(res=pth)

			source = fi.FITS(pth)
			matched_hoopoe.res = source["i3s_omhoopoe"].read()
			matched_hoopoe.truth = source["truth_omhoopoe"].read()
			matched_waxwing.res = source["i3s_omwaxwing"].read()
			matched_waxwing.truth = source["truth_omwaxwing"].read()

			matched_hoopoe1 = s.shapecat()
			matched_waxwing1 = s.shapecat()
			matched_hoopoe2 = s.shapecat()


			cut1=((matched_hoopoe.res["snr"] > 12) & (matched_hoopoe.res["snr"] < 200) & (matched_hoopoe.res["mean_rgpp_rp"] > 1.13) & (matched_hoopoe.res["mean_rgpp_rp"] < 3.0) & (matched_hoopoe.res["info_flag"]==0)) 
			hoopoe.res=matched_hoopoe.res[cut1]
			hoopoe.truth=matched_hoopoe.truth[cut1]

			cut2=((matched_waxwing.res["snr"] > 12) & (matched_waxwing.res["snr"] < 200) & (matched_waxwing.res["mean_rgpp_rp"] > 1.13) & (matched_waxwing.res["mean_rgpp_rp"] < 3.0) & (matched_waxwing.res["info_flag"]==0)) 
			waxwing.res=matched_waxwing.res[cut2]
			waxwing.truth=matched_waxwing.truth[cut2]

			matched_waxwing1.res=matched_waxwing.res[cut2 & cut1]
			matched_waxwing1.truth=matched_waxwing.truth[cut2 & cut1]

			matched_hoopoe1.res=matched_hoopoe.res[cut2 & cut1]
			matched_hoopoe1.truth=matched_hoopoe.truth[cut2 & cut1]

			hoopoe_diluted=s.shapecat()

			hoopoe_diluted.res = matched_hoopoe.res[cut1 & cut2]
			hoopoe_diluted.truth = matched_hoopoe.truth[cut1 & cut2]
			max_frac_bin_shift=0.22
			select_hoopoe_not_waxwing = ((abs(matched_waxwing.res["mean_rgpp_rp"]-matched_hoopoe.res["mean_rgpp_rp"])/matched_hoopoe.res["mean_rgpp_rp"] >max_frac_bin_shift) | (abs(np.log10(matched_waxwing.res["snr"])-np.log10(matched_hoopoe.res["snr"]))/np.log10(matched_hoopoe.res["snr"]) >max_frac_bin_shift) ) & cut1
			#select_hoopoe_not_waxwing = ((abs(np.log10(matched_waxwing.res["snr"])-np.log10(matched_hoopoe.res["snr"]))/np.log10(matched_waxwing.res["snr"]) >0.18) ) & cut1
			#select_hoopoe_not_waxwing = ((abs(matched_waxwing.res["mean_rgpp_rp"]-matched_hoopoe.res["mean_rgpp_rp"])/matched_waxwing.res["mean_rgpp_rp"] >0.24) ) & cut1
			import copy
			hoopoe_not_waxwing = copy.deepcopy(matched_hoopoe.res[select_hoopoe_not_waxwing])
			hoopoe_not_waxwing["e1"] = np.random.normal(0,1, hoopoe_not_waxwing["e1"].size)*hoopoe_not_waxwing["e1"].std()
			hoopoe_not_waxwing["e2"] = np.random.normal(0,1, hoopoe_not_waxwing["e2"].size)*hoopoe_not_waxwing["e2"].std()
			hoopoe_diluted.res = np.hstack((hoopoe_diluted.res,hoopoe_not_waxwing))

			hoopoe_not_waxwing = copy.deepcopy(matched_hoopoe.truth[select_hoopoe_not_waxwing])
			#hoopoe_not_waxwing["true_g1"] = np.random.rand(hoopoe_not_waxwing["true_g1"].size)*0.16-0.08
			#hoopoe_not_waxwing["true_g2"] = np.random.rand(hoopoe_not_waxwing["true_g2"].size)*0.16-0.08
			hoopoe_diluted.truth = np.hstack((hoopoe_diluted.truth,hoopoe_not_waxwing))


			b1, b2 = np.log10(12), np.log10(200)
			vec_matched1_h_diluted = hoopoe_diluted.bias_vs_obs("snr", "m", 9, binning=np.logspace(b1,b2,10),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--")
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-snr-hoopoe-matched-pergal_diluted_randomise_eobs",np.vstack((vec_matched1_h_diluted[0], vec_matched1_h_diluted[1]["m"], vec_matched1_h_diluted[2]["m"])).T)
			vec_unmatched_h = hoopoe.bias_vs_obs("snr", "m", 9, binning=np.logspace(b1,b2,10),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--")
			vec_unmatched_nf = waxwing.bias_vs_obs("snr", "m", 9, binning=np.logspace(b1,b2,10),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--")
			vec_matched_h = matched_hoopoe1.bias_vs_obs("snr", "m", 9, binning=np.logspace(b1,b2,10),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--", xdata=matched_waxwing1.res["snr"])
			vec_matched_nf = matched_waxwing1.bias_vs_obs("snr", "m", 9, binning=np.logspace(b1,b2,10),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--", xdata=matched_waxwing1.res["snr"])
			vec_matched1_h = matched_hoopoe1.bias_vs_obs("snr", "m", 9, binning=np.logspace(b1,b2,10),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--")
			vec_matched2_nf = matched_waxwing1.bias_vs_obs("snr", "m", 9, binning=np.logspace(b1,b2,10),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--")

			if alpha:
				cvec_unmatched_h=hoopoe.bias_vs_obs("snr", "c", 9, binning=np.logspace(b1,b2,10),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--")
				cvec_unmatched_nf=waxwing.bias_vs_obs("snr", "c", 9, binning=np.logspace(b1,b2,10),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--")
				cvec_matched_h=matched_hoopoe1.bias_vs_obs("snr", "c", 9, binning=np.logspace(b1,b2,10),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--", xdata=matched_waxwing1.res["snr"])
				cvec_matched_nf=matched_waxwing1.bias_vs_obs("snr", "c", 9, binning=np.logspace(b1,b2,10),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--", xdata=matched_waxwing1.res["snr"])
				cvec_matched2_h=matched_hoopoe1.bias_vs_obs("snr", "c", 9, binning=np.logspace(b1,b2,10),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--")
				cvec_matched2_nf=matched_waxwing1.bias_vs_obs("snr", "c", 9, binning=np.logspace(b1,b2,10),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--")

				nbins_alpha=6
				avec_unmatched_h=hoopoe.alpha_vs_obs("snr", "alpha", nbins_alpha, binning=np.logspace(b1,b2,nbins_alpha+1),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--")
				avec_unmatched_nf=waxwing.alpha_vs_obs("snr", "alpha", nbins_alpha, binning=np.logspace(b1,b2,nbins_alpha+1),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--")
				avec_matched_h=matched_hoopoe1.alpha_vs_obs("snr", "alpha", nbins_alpha, binning=np.logspace(b1,b2,nbins_alpha+1),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--", xdata=matched_waxwing1.res["snr"])
				avec_matched_nf=matched_waxwing1.alpha_vs_obs("snr", "alpha", nbins_alpha, binning=np.logspace(b1,b2,nbins_alpha+1),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--", xdata=matched_waxwing1.res["snr"])
				avec_matched2_h=matched_hoopoe1.alpha_vs_obs("snr", "alpha", nbins_alpha, binning=np.logspace(b1,b2,nbins_alpha+1),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--")
				avec_matched2_nf=matched_waxwing1.alpha_vs_obs("snr", "alpha", nbins_alpha, binning=np.logspace(b1,b2,nbins_alpha+1),error_type="std", logscale=True, colour="purple", return_vals=True, ls="--")

			import pdb ; pdb.set_trace()

			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-snr-hoopoe-unmatched",np.vstack((vec_unmatched_h[0], vec_unmatched_h[1]["m"], vec_unmatched_h[2]["m"])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-snr-hoopoe-matched",np.vstack((vec_matched_h[0], vec_matched_h[1]["m"], vec_matched_h[2]["m"])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-snr-waxwing-matched",np.vstack((vec_matched_nf[0], vec_matched_nf[1]["m"], vec_matched_nf[2]["m"])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-snr-waxwing-unmatched",np.vstack((vec_unmatched_nf[0], vec_unmatched_nf[1]["m"], vec_unmatched_nf[2]["m"])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-snr-waxwing-matched-pergal",np.vstack((vec_matched2_nf[0], vec_matched2_nf[1]["m"], vec_matched2_nf[2]["m"])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-snr-hoopoe-matched-pergal",np.vstack((vec_matched1_h[0], vec_matched1_h[1]["m"], vec_matched1_h[2]["m"])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-snr-hoopoe-matched-pergal_diluted_randomise_eobs",np.vstack((vec_matched1_h_diluted[0], vec_matched1_h_diluted[1]["m"], vec_matched1_h_diluted[2]["m"])).T)

			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/a-vs-snr-hoopoe-unmatched",np.vstack((avec_unmatched_h[0], avec_unmatched_h[1], avec_unmatched_h[2])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/a-vs-snr-hoopoe-matched",np.vstack((avec_matched_h[0], avec_matched_h[1], avec_matched_h[2])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/a-vs-snr-waxwing-matched",np.vstack((avec_matched_nf[0], avec_matched_nf[1], avec_matched_nf[2])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/a-vs-snr-waxwing-unmatched",np.vstack((avec_unmatched_nf[0], avec_unmatched_nf[1], avec_unmatched_nf[2])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/a-vs-snr-waxwing-matched-pergal",np.vstack((avec_matched2_nf[0], avec_matched2_nf[1], avec_matched2_nf[2])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/a-vs-snr-hoopoe-matched-pergal",np.vstack((avec_matched2_h[0], avec_matched2_h[1], avec_matched2_h[2])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/a-vs-snr-hoopoe-unmatched",np.vstack((avec_unmatched_h[0], avec_unmatched_h[1], avec_unmatched_h[2])).T)

		if generate_data_r:
			pth="/share/des/disc8/cambridge/bord-fits/combined_cats/matched_hoopoe-waxwing_nocuts.fits"
			pth="/share/des/disc8/y1-unblinded-combined_cats/hoopoe-waxwing-overlap-unblinded.fits"
			print pth
			matched_hoopoe =s.shapecat(res=pth)
			matched_waxwing =s.shapecat(res=pth)
			hoopoe =s.shapecat(res=pth)
			waxwing =s.shapecat(res=pth)

			source = fi.FITS(pth)
			matched_hoopoe.res = source["i3s_omhoopoe"].read()
			matched_hoopoe.truth = source["truth_omhoopoe"].read()
			matched_waxwing.res = source["i3s_omwaxwing"].read()
			matched_waxwing.truth = source["truth_omwaxwing"].read()

			matched_hoopoe1 = s.shapecat()
			matched_waxwing1 = s.shapecat()
			matched_hoopoe2 = s.shapecat()


			cut1=((matched_hoopoe.res["snr"] > 12) & (matched_hoopoe.res["snr"] < 200) & (matched_hoopoe.res["mean_rgpp_rp"] > 1.13) & (matched_hoopoe.res["mean_rgpp_rp"] < 1.8) & (matched_hoopoe.res["info_flag"]==0)) 
			hoopoe.res=matched_hoopoe.res[cut1]
			hoopoe.truth=matched_hoopoe.truth[cut1]

			cut2=((matched_waxwing.res["snr"] > 12) & (matched_waxwing.res["snr"] < 200) & (matched_waxwing.res["mean_rgpp_rp"] > 1.13) & (matched_waxwing.res["mean_rgpp_rp"] < 1.8) & (matched_waxwing.res["info_flag"]==0)) 
			waxwing.res=matched_waxwing.res[cut2]
			waxwing.truth=matched_waxwing.truth[cut2]

			matched_waxwing1.res=matched_waxwing.res[cut2 & cut1]
			matched_waxwing1.truth=matched_waxwing.truth[cut2 & cut1]

			matched_hoopoe1.res=matched_hoopoe.res[cut2 & cut1]
			matched_hoopoe1.truth=matched_hoopoe.truth[cut2 & cut1]

			hoopoe_diluted = s.shapecat()

			hoopoe_diluted.res = matched_hoopoe.res[cut1 & cut2]
			hoopoe_diluted.truth = matched_hoopoe.truth[cut1 & cut2]
			max_frac_bin_shift=0.22
			select_hoopoe_not_waxwing = ((abs(matched_waxwing.res["mean_rgpp_rp"]-matched_hoopoe.res["mean_rgpp_rp"])/matched_hoopoe.res["mean_rgpp_rp"] >max_frac_bin_shift) | (abs(np.log10(matched_waxwing.res["snr"])-np.log10(matched_hoopoe.res["snr"]))/np.log10(matched_hoopoe.res["snr"]) >max_frac_bin_shift) ) & cut1
			#select_hoopoe_not_waxwing = ((abs(np.log10(matched_waxwing.res["snr"])-np.log10(matched_hoopoe.res["snr"]))/np.log10(matched_waxwing.res["snr"]) >0.18) ) & cut1
			#select_hoopoe_not_waxwing = ((abs(matched_waxwing.res["mean_rgpp_rp"]-matched_hoopoe.res["mean_rgpp_rp"])/matched_waxwing.res["mean_rgpp_rp"] >0.24) ) & cut1
			import copy
			hoopoe_not_waxwing = copy.deepcopy(matched_hoopoe.res[select_hoopoe_not_waxwing])
			hoopoe_not_waxwing["e1"] = np.random.normal(0,1, hoopoe_not_waxwing["e1"].size)*hoopoe_not_waxwing["e1"].std()
			hoopoe_not_waxwing["e2"] = np.random.normal(0,1, hoopoe_not_waxwing["e2"].size)*hoopoe_not_waxwing["e2"].std()
			hoopoe_diluted.res = np.hstack((hoopoe_diluted.res,hoopoe_not_waxwing))

			hoopoe_not_waxwing = copy.deepcopy(matched_hoopoe.truth[select_hoopoe_not_waxwing])
			#hoopoe_not_waxwing["true_g1"] = np.random.rand(hoopoe_not_waxwing["true_g1"].size)*0.16-0.08
			#hoopoe_not_waxwing["true_g2"] = np.random.rand(hoopoe_not_waxwing["true_g2"].size)*0.16-0.08
			hoopoe_diluted.truth = np.hstack((hoopoe_diluted.truth,hoopoe_not_waxwing))



			#dilution_mask = np.invert(cut2) & cut1
			#hoopoe_diluted = s.shapecat()
			#hoopoe_diluted.res = matched_hoopoe.res[cut1]
			#hoopoe_diluted.truth = matched_hoopoe.truth
			#nrandomise = hoopoe_diluted.truth["true_g1"][dilution_mask].size
			#hoopoe_diluted.truth["true_g1"][dilution_mask] = np.random.rand(nrandomise)*0.16-0.08
			#hoopoe_diluted.truth["true_g2"][dilution_mask] = np.random.rand(nrandomise)*0.16-0.08
			#remove_again = dilution_mask & (np.random.rand(dilution_mask.size)>0.3333)
			#print nrandomise, remove_again.size
			#hoopoe_diluted.truth["true_g1"][remove_again]=1e8 
			#hoopoe_diluted.truth["true_g2"][remove_again]=1e8 
			#print hoopoe_diluted.truth["true_g1"][remove_again]
			#hoopoe_diluted.truth = hoopoe_diluted.truth[cut1]

			b1, b2 = 1.13, 1.7
			nrbin=9
			bin_edges=di.find_bin_edges(matched_waxwing1.res["mean_rgpp_rp"], nrbin)
			vec_hoopoe_diluted = hoopoe_diluted.bias_vs_obs("mean_rgpp_rp", "m", nrbin, binning=bin_edges,error_type="std", logscale=True, colour="purple", return_vals=True, ls="--")
			import pdb ; pdb.set_trace()
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-rgpp-hoopoe-matched-pergal_diluted_randomised_eobs",np.vstack((vec_hoopoe_diluted[0], vec_hoopoe_diluted[1]["m"], vec_hoopoe_diluted[2]["m"])).T)
			vec_hoopoe_owncuts = hoopoe.bias_vs_obs("mean_rgpp_rp", "m", nrbin, binning=bin_edges,error_type="std", logscale=False, colour="purple", return_vals=True, ls="--")
			vec_waxwing_owncuts = waxwing.bias_vs_obs("mean_rgpp_rp", "m", nrbin, binning=bin_edges,error_type="std", logscale=False, colour="purple", return_vals=True, ls="--")
			vec_hoopoe_bins_matched = matched_hoopoe1.bias_vs_obs("mean_rgpp_rp", "m", nrbin, binning=bin_edges,error_type="std", logscale=False, colour="purple", return_vals=True, ls="--", xdata=matched_waxwing1.res["mean_rgpp_rp"])
			vec_waxwing_bins_matched = matched_waxwing1.bias_vs_obs("mean_rgpp_rp", "m", nrbin, binning=bin_edges,error_type="std", logscale=False, colour="purple", return_vals=True, ls="--", xdata=matched_waxwing1.res["mean_rgpp_rp"])
			vec_hoopoe_objects_matched = matched_hoopoe1.bias_vs_obs("mean_rgpp_rp", "m", nrbin, binning=bin_edges,error_type="std", logscale=False, colour="purple", return_vals=True, ls="--")
			vec_matched2_nf = matched_waxwing1.bias_vs_obs("mean_rgpp_rp", "m", nrbin, binning=bin_edges,error_type="std", logscale=False, colour="purple", return_vals=True, ls="--")


			nbins_alpha=6
			if alpha:
				avec_unmatched_h=hoopoe.alpha_vs_obs("mean_rgpp_rp", "alpha", nbins_alpha, binning="equal_number",error_type="std", logscale=False, colour="purple", return_vals=True, ls="--")
				avec_unmatched_nf=waxwing.alpha_vs_obs("mean_rgpp_rp", "alpha", nbins_alpha, binning="equal_number",error_type="std", logscale=False, colour="purple", return_vals=True, ls="--")
				avec_matched_h=matched_hoopoe1.alpha_vs_obs("mean_rgpp_rp", "alpha", nbins_alpha, binning="equal_number",error_type="std", logscale=False, colour="purple", return_vals=True, ls="--", xdata=matched_waxwing1.res["mean_rgpp_rp"])
				avec_matched_nf=matched_waxwing1.alpha_vs_obs("mean_rgpp_rp", "alpha", nbins_alpha, binning="equal_number",error_type="std", logscale=False, colour="purple", return_vals=True, ls="--", xdata=matched_waxwing1.res["mean_rgpp_rp"])
				avec_matched2_h=matched_hoopoe1.alpha_vs_obs("mean_rgpp_rp", "alpha", nbins_alpha, binning="equal_number",error_type="std", logscale=False, colour="purple", return_vals=True, ls="--")
				avec_matched2_nf=matched_waxwing1.alpha_vs_obs("mean_rgpp_rp", "alpha", nbins_alpha, binning="equal_number",error_type="std", logscale=False, colour="purple", return_vals=True, ls="--")

				np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/a-vs-rgpp-hoopoe-unmatched",np.vstack((avec_unmatched_h[0], avec_unmatched_h[1], avec_unmatched_h[2])).T)
				np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/a-vs-rgpp-hoopoe-matched",np.vstack((avec_matched_h[0], avec_matched_h[1], avec_matched_h[2])).T)
				np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/a-vs-rgpp-waxwing-matched",np.vstack((avec_matched_nf[0], avec_matched_nf[1], avec_matched_nf[2])).T)
				np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/a-vs-rgpp-waxwing-unmatched",np.vstack((avec_unmatched_nf[0], avec_unmatched_nf[1], avec_unmatched_nf[2])).T)
				np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/a-vs-rgpp-waxwing-matched-pergal",np.vstack((avec_matched2_nf[0], avec_matched2_nf[1], avec_matched2_nf[2])).T)
				np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/a-vs-rgpp-hoopoe-matched-pergal",np.vstack((avec_matched2_h[0], avec_matched2_h[1], avec_matched2_h[2])).T)
				np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/a-vs-rgpp-hoopoe-unmatched",np.vstack((avec_unmatched_h[0], avec_unmatched_h[1], avec_unmatched_h[2])).T)

			import pdb ; pdb.set_trace()
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-rgpp-hoopoe-matched-pergal_diluted",np.vstack((vec_hoopoe_diluted[0], vec_hoopoe_diluted[1]["m"], vec_hoopoe_diluted[2]["m"])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-rgpp-hoopoe-unmatched",np.vstack((vec_hoopoe_owncuts[0], vec_hoopoe_owncuts[1]["m"], vec_hoopoe_owncuts[2]["m"])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-rgpp-hoopoe-matched",np.vstack((vec_hoopoe_bins_matched[0], vec_hoopoe_bins_matched[1]["m"], vec_hoopoe_bins_matched[2]["m"])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-rgpp-waxwing-matched",np.vstack((vec_waxwing_bins_matched[0], vec_waxwing_bins_matched[1]["m"], vec_waxwing_bins_matched[2]["m"])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-rgpp-waxwing-unmatched",np.vstack((vec_waxwing_owncuts[0], vec_waxwing_owncuts[1]["m"], vec_waxwing_owncuts[2]["m"])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-rgpp-waxwing-matched-pergal",np.vstack((vec_matched2_nf[0], vec_matched2_nf[1]["m"], vec_matched2_nf[2]["m"])).T)
			np.savetxt("/home/samuroff/hoopoe_paper/plots/data/fig4/m-vs-rgpp-hoopoe-matched-pergal",np.vstack((vec_hoopoe_objects_matched[0], vec_hoopoe_objects_matched[1]["m"], vec_hoopoe_objects_matched[2]["m"])).T)

	def summary_plot(self):

		plt.close()
		plt.switch_backend("pdf")
		cases = ["hf", "wf", "m_n1f1_owncuts", "m_n0f0_owncuts", "m_n0f0_hoopoecuts", "m_n1f1_owncuts+dgn", "m_n1f1_bothcuts+dgn", "m_n0f0_bothcuts", "m2_n1f0_owncuts", "m2_n1f1_owncuts", "m2_n1f1_f0cuts"]

#		lookup={
#		"hf":[-0.15134869469909246, 0.001831402995716225, "$hoopoe$, own cuts"],
#		"wf":[-0.099733818389343892, 0.0049929106165088272, "$waxwing$, own cuts"],
#		"m_n1f1_owncuts":[-0.14605422357680206, 0.0050273319928490897, "Matched $hoopoe$, own cuts"],
#		"m_n1f1_owncuts+dgn":[-0.11333853344766834, 0.0063829154847302736, r"Matched $hoopoe$, own cuts \& $d_{gn}>20$ pix"],
#		"m_n1f1_bothcuts+dgn":[-0.10544, 0.005697, r"Matched $hoopoe$, both cuts \& $d_{gn}>20$ pix"],
#		"m_n0f0_owncuts":[-0.099374117849417298, 0.0052202293676924597, "Matched $waxwing$, own cuts"],
#		"m_n0f0_hoopoecuts":[-0.089707741757402717, 0.0051494746527769777, "Matched $waxwing$, $hoopoe$ cuts"],
#		"m_n0f0_bothcuts":[-0.103985, 0.0058428, "Matched $waxwing$, both cuts \& $d_{gn}>20$ pix"],
#		"m2_n1f0_owncuts":[-0.13806,0.005216, "No sub-detection, own cuts" ],
#		"m2_n1f1_owncuts":[-0.1539,0.005842, "Matched $hoopoe$, own cuts" ],
#		"m2_n1f1_f0cuts":[-0.1341,0.00753, "$hoopoe$, cuts from no sub-detection" ]
#		}

		lookup={
		"hf":[-0.11902291590209725, 0.0017383162244452117, "$hoopoe$, own cuts"],
		"wf":[-0.064, 0.006, "$waxwing$, own cuts"],
		"m_n1f1_owncuts":[-0.113, 0.007, "Matched $hoopoe$, own cuts"],
		"m_n1f1_owncuts+dgn":[-0.078, 0.008, r"Matched $hoopoe$, own cuts & $d_{gn}>20$ pix"],
		"m_n1f1_bothcuts+dgn":[-0.070, 0.007, r"Matched $hoopoe$, both cuts & $d_{gn}>20$ pix"],
		"m_n0f0_owncuts":[-0.064, 0.006, "Matched $waxwing$, own cuts"],
		"m_n0f0_hoopoecuts":[-0.055, 0.007, r"Matched $waxwing$, $hoopoe$ cuts"],
		"m_n0f0_bothcuts":[-0.068, 0.007, r"Matched $waxwing$, both cuts & $d_{gn}>20$ pix"],
		"m2_n1f0_owncuts":[-0.104, 0.008, "No sub-detection, own cuts" ],
		"m2_n1f1_owncuts":[-0.120, 0.006, "Matched $hoopoe$, own cuts" ],
		"m2_n1f1_f0cuts":[-0.101, 0.009, "$hoopoe$, cuts from no sub-detection" ]
		}

		colours=["purple", "steelblue", "forestgreen", "pink", "darkred", "orange", "k", "darkviolet"]
		fontcolours=["darkviolet"]*2 + ["forestgreen"]*6 + ["red"]*3

		fig = plt.figure(1)

		plt.ylim(100*(1+len(cases)),0.0)
		plt.xlim(-0.175,-0.045)
		ax = plt.gca()
		ax.yaxis.set_visible(False)

		y = np.arange(100,100*(1+len(cases)),100)
		for i, (y0,name) in enumerate(zip(y,cases)):
			line = lookup[name]
			print i, line
			if i==0:
				pt="*" ; fc="none" ; sz=9 ; ec="k"
			else:
				pt="x" ; fc=fontcolours[i] ; sz=6 ; ec=fontcolours[i]
			plt.errorbar(line[0], y0, xerr=line[1], color=fontcolours[i], mec=ec, fmt=pt, markerfacecolor=fc, elinewidth=1.5, markersize=sz)
			plt.text(-0.174, y0, line[2], fontsize=8.4, color=fontcolours[i])
			if i==0:
				plt.axvline(line[0], color=colours[i])
				plt.axvspan(line[0]-line[1], line[0]+line[1], color="purple", alpha=0.2)

		plt.xlabel("Multiplicative Bias $m \equiv (m_1+m_2)/2$")
		plt.subplots_adjust(top=0.98,bottom=0.138)

		plt.savefig("/home/samuroff/shear_pipeline/plot_dump/hoopoe/m_summary_plot-unblinded.pdf")


#import tools.diagnostics as di
#import tools.plots as pl
#import fitsio as fi
#import pylab as plt
#import os, yaml, argparse
#import tools.shapes as s
#import tools.nbc as cal
#import tools.arrays as arr
#from tools.im3shape import calibrate_all_tomographic as ct
#from tools.im3shape import calibrate_all as ca
#import copy
#from tools.im3shape.calibration import calibrate_y1 as cy1
#config = yaml.load(open("/home/samuroff/local/python/lib/python2.7/site-packages/tools/im3shape/calibration/config/fiducial-y1-unblinded.yaml"))
#hoopoe, weights, y1v2 = cy1.setup(True, False, config)
#
#
## exposure map
#
#plt.style.use("y1a1")
#plt.switch_backend("pdf")
#col=hoopoe.map_depth(thin=3)
#sky_map(hoopoe.res["ra"][::3],hoopoe.res["dec"][::3],colour=col["n_exposure"])
#
#
## Toy model shear resiudal
#
#dv=np.genfromtxt("/home/samuroff/shear_residual-toymodel-median_vals-dgn6.txt",names=True)
#
#dv0=np.genfromtxt("/home/samuroff/shear_residual-toymodel-median_vals-dgn6-noneigh.txt",names=True)
#
#plt.close()
#import matplotlib
#matplotlib.rcParams['font.family']='serif'
#matplotlib.rcParams['font.size']=18
#matplotlib.rcParams['legend.fontsize']=18
#
#
#plt.close()
#fig, ax1 = plt.subplots()
##ax2 = fig.add_axes([0.27, 0.205, 0.265, 0.265])
##ax3 = fig.add_axes([0.72, 0.7, 0.25, 0.25])
#ax1.plot(dv["g"], dv["dg"],"o", color="purple",mec="purple")
#ax1.plot(dv0["g"], dv0["dg"],":", color="purple")
#ax1.axhline(0, color="k")
#ax1.axvline(0, color="k")
#ax1.set_ylim(-0.09,0.09)
#ax1.set_xlim(-0.31,0.31)
#ax1.set_yticks([-0.09,-0.06,-0.03,0,0.03,0.06,0.09])
#ax1.plot([-0.6,0.6], np.array([-0.6,0.6])*-0.45, lw=1, color="midnightblue", alpha=0.8)
#ax1.plot([-0.6,0.6], np.array([-0.6,0.6])*-0.5, lw=1, color="midnightblue", alpha=0.8)
#ax1.plot([-0.6,0.6], np.array([-0.6,0.6])*-0.4, lw=1, color="midnightblue", alpha=0.8)
#ax1.set_xlabel(r"Input Shear $g^{tr}$", fontsize=22)
#ax1.set_ylabel(r"Shear Residual $\tilde{g}-g^{tr}$", fontsize=22)

#ax2.plot(dv["g"], dv["dg"],".", color="purple")
#ax2.plot(dv0["g"], dv0["dg"],":", color="purple")
#ax2.axhline(0, color="k")
#ax2.axvline(0, color="k")
#ax2.set_ylim(-0.04,0.04)
#ax2.set_xlim(-0.08,0.08)
#ax2.set_yticks([-0.04,-0.02,0.00,0.02,0.04])
#ax2.set_xticks([-0.08, -0.04,0.00,0.04,0.08])
#ax2.plot([-0.6,0.6], np.array([-0.6,0.6])*-0.45, lw=1, color="k", alpha=0.8)
#ax2.plot([-0.6,0.6], np.array([-0.6,0.6])*-0.5, lw=1, color="k", alpha=0.8)
#ax2.plot([-0.6,0.6], np.array([-0.6,0.6])*-0.4, lw=1, color="k", alpha=0.8)
#ax2.tick_params(labelsize=13)

plt.subplots_adjust(left=0.2, right=0.98, top=0.97,bottom=0.15)

plt.savefig("/home/samuroff/g-vs-dg-toy_model_medians.pdf")
plt.close()

#
#
## cartoon ellipses
#
#
#from tools.im3shape import basic as i3s
#import tools.shapes as s
#m=s.meds_wrapper("/share/des/disc8/cambridge/meds/DES2111+0043-r-sim-ohioA6-meds-y1a1-beta.fits.fz")
#g_r0={}
#r=np.linspace(0,35,60)
#reload(pl)
#
#
#params={
#    "Rc":2.1,
#    "fc":945,
#    "fn":475,
#    "Rn":1.47,
#    "Rp":0.8}
#
#plt.close() 
#plt.style.use("y1a1")
#plt.switch_backend("pdf")
#plt.gca(aspect="equal", adjustable="box")
#gal0,psf0=i3s.setup_simple(boxsize=32*3,shear=(-0.4,0.0), psf_size=params["Rp"],  size=params["Rc"]*1.5*3, neighbour_ellipticity=(0.0,0.0), neighbour_flux=params["fn"], flux=params["fc"]*3, neighbour_size=0.5, neighbour=[np.inf,np.inf], opt=m.options)
#dt=np.pi/8.
#np.cos(dt)
#dgn=6*3
#dx,dy=3*np.cos(dt),3*np.sin(dt)
#xy=np.meshgrid(np.arange(-32*3/2,32*3/2,1),np.arange(-32*3/2,32*3/2,1) )
#CS = plt.contour(xy[0]-0.5,xy[1]-0.5,gal0,levels=np.linspace(0.115,13,10), colors="purple")
#plt.scatter([dgn],[0], marker="x", s=90, color="steelblue")
#plt.scatter([0],[dgn], marker="x", s=90, color="steelblue")
#plt.scatter([dgn/np.sqrt(2)],[dgn/np.sqrt(2)], marker="x", s=90, color="steelblue")
#plt.arrow(0,dgn,-3,0, color="steelblue", head_width=0.65, head_length=0.75, lw=1.5)
#plt.arrow(dgn/np.sqrt(2),dgn/np.sqrt(2), -3/2,3/2, color="steelblue", head_width=0.65, head_length=0.75, lw=1.5)
#plt.arrow(dgn,0,0,3, color="steelblue", head_width=0.65, head_length=0.75, lw=1.5)
#plt.text(dgn+1,0, "A", fontsize=26)
#plt.text(dgn/np.sqrt(2)+1, dgn/np.sqrt(2)+1, "B", fontsize=26)
#plt.text(0, dgn+1, "C", fontsize=26)
#plt.scatter([0],[0],marker="+", color="k", s=200,lw=2)
#plt.xlim(-31,30)
#plt.ylim(-31,30)
#plt.xticks(visible=False)
#plt.yticks(visible=False)
#
#plt.savefig("tmp12.pdf")
#
#
#
## toy model et vs dgn
#
#plt.style.use("y1a1")
#plt.switch_backend("pdf")
#plt.close()
#import matplotlib
#matplotlib.rcParams['font.family']='serif'
#matplotlib.rcParams['legend.fontsize']=22
#matplotlib.rcParams['font.size']=18
#
#r,g1,g2 = np.loadtxt("/home/samuroff/etangent-vs-dgn-toy_model_medians-round_neighbour.txt").T
#r,g1n,g2n = np.loadtxt("/home/samuroff/etangent-vs-dgn-toy_model_medians-sheared_neighbour.txt").T
#
#
#plt.plot(r, -1*np.array(g1), color="purple", label=r"$\tilde{g}_1$ (round neighbour)", lw=2.5, ls="-")
#plt.plot(r, -1*np.array(g2), color="purple", label=r"$\tilde{g}_2$ (round neighbour)", lw=2.5, ls="-.")
#plt.plot(r, -1*np.array(g1n), color="k", label=r"$\tilde{g}_1$ ($g^{tr}_{2,n}=0.1$)", lw=2.5, ls=":")
#plt.plot(r, -1*np.array(g2n), color="k", label=r"$\tilde{g}_2$ ($g^{tr}_{2,n}=0.1$)", lw=2.5, ls="--")
#plt.axhline(0, color="k")
#plt.xlim(0.,32)
#plt.ylim(-0.45,0.04)
#plt.legend(loc="lower right")
#plt.xlabel("Neighbour Distance $d_{gn}$ / pixels", fontsize=22)
#plt.ylabel(r"Tangential or Cross Shear", fontsize=22)
##$- \tilde{g}_i (d_{gn}| g^{tr}=0, \theta=0) $")
#plt.subplots_adjust(left=0.16, right=0.98, bottom=0.15, top=0.98)
#plt.savefig("/home/samuroff/dgn-vs-etangent-toy_model_medians.pdf")
#
#
#
#def line(x,a,b):
#    return np.log10(x)*a+b
#
#def curve(x,a,b):
#    return b/(x**a)
#
## gm and mm correlation functions
#plt.close()
#plt.style.use("y1a1")
#plt.switch_backend("pdf")
#import matplotlib
#matplotlib.rcParams["legend.fontsize"]=22
#
#m0 =  0.15126887361063673
#xmm, mm = np.loadtxt("/home/samuroff/hoopoe_paper/paper/data/2pt_m/2pt-mm-binslop_0.100.txt").T
#xgm, gm = np.loadtxt("/home/samuroff/hoopoe_paper/paper/data/2pt_m/2pt-gm-binslop_0.100.txt").T
#
#bf_mm,cov_mm = sp.optimize.curve_fit(curve, xmm, mm-m0*m0, p0=[1,0.02])
#bf_gm,cov_gm = sp.optimize.curve_fit(line, xgm, gm)
#
#plt.close()
#plt.plot(xmm, (mm-m0*m0)*1e3, "o", color="purple", mec="none", label=r"$\left < m m \right > - \bar{m}^2$")
#plt.plot(xgm, gm*1e3, "^", color="steelblue", mec="none", label=r"$\left < \delta_g m \right >$")
#plt.plot(np.logspace(0,3,400), 1e3*curve(np.logspace(0,3,400), bf_mm[0], bf_mm[1]), ":", color="purple", lw=2.5)
#plt.plot(np.logspace(0,3,400), 1e3*line(np.logspace(0,3,400), bf_gm[0], bf_gm[1]), ":", color="steelblue", lw=2.5)
#plt.axvline(0.75/6*60*np.sqrt(2),color="k", ls="--")
#plt.axhline(0, color="k")
#plt.axvspan(1,3.6, color="blue", alpha=0.1)
#plt.axvspan(250,500, color="blue", alpha=0.1)
#plt.xlim(2,250)
#plt.legend(numpoints=3)
#plt.xlabel(r"Angular Scale $\theta$ / arcseconds",fontsize=22)
#plt.ylabel(r"$\xi_{a b}(\theta)$ / $\times10^{-3}$", fontsize=22)
#plt.ylim(-4.5,4.5)
#plt.xscale("log")
#plt.subplots_adjust(bottom=0.13)
#plt.savefig("/home/samuroff/shear_pipeline/plot_dump/hoopoe/xi-mm-gm-6x6-fit-v2.pdf")
#
#plt.close()
#plt.plot(xgg, gg*1e3, "<", color="deeppink", mec="none", label=r"$\left < \delta_g \delta_g \right >$")
#plt.plot(np.logspace(0,3,400), 1e3*curve(np.logspace(0,3,400), bf_gg[0], bf_gg[1]), ":", color="deeppink")
#plt.axvline(0.75/6*60*np.sqrt(2),color="k", ls="--")
#plt.axhline(0, color="k")
#plt.axvspan(1,3.6, color="blue", alpha=0.1)
#plt.axvspan(250,500, color="blue", alpha=0.1)
#plt.xlim(2,250)
#plt.legend()
#plt.xlabel(r"Angular Scale $\theta$ / arcseconds",fontsize=22)
#plt.ylabel(r"$\xi_{a b}(\theta)$ / $\times10^{-3}$", fontsize=22)
#plt.xscale("log")
#plt.subplots_adjust(bottom=0.13)
#plt.savefig("/home/samuroff/shear_pipeline/plot_dump/hoopoe/xi-delta_g-delta_g-6x6-fit.pdf")





