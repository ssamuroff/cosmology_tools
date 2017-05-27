import numpy as np
import tools.nbc as cal
import tools.shapes as s
import tools.diagnostics as di
import tools.plots as pl
import tools.arrays as arr
import fitsio as fi
import pylab as plt
import os, yaml, argparse, glob, gc

plt.switch_backend("pdf")

def setup(load_sim, load_data, config, verbose=True):
	im3shape_columns = ["e1", "e2", "mean_hsm_psf_e1_sky", "mean_hsm_psf_e2_sky", "mean_hsm_psf_sigma", "snr", "mean_rgpp_rp","mean_mask_fraction", "radius", "coadd_objects_id", "is_bulge", "bulge_flux", "disc_flux", "info_flag", "mag_auto_r"]
	truth_columns = ['DES_id', 'cosmos_ident', 'cosmos_photoz', 'sextractor_pixel_offset', 'true_g1', 'true_g2', 'intrinsic_e1', 'intrinsic_e2', 'ra', 'dec', 'hlr', 'mag', 'flux']

	y1v2 = None
	hoopoe = None

	# Load the y1 data
	if load_data:
		if verbose:
			print "Loading data %s"%config["input"]["i3s"]
		y1v2 = s.shapecat(res=config["input"]["i3s"])
		y1v2.load(truth=False, prune=True, cols=[im3shape_columns,truth_columns])
		y1v2.res = y1v2.res[y1v2.res["info_flag"]==0] # This should always be true, but just in case...
		sel = ((y1v2.res["snr"] > 12) & (y1v2.res["snr"] < 200) & (y1v2.res["mean_rgpp_rp"] > 1.13) & (y1v2.res["mean_rgpp_rp"] < 3.0))
		y1v2.res=y1v2.res[sel]
	else:
		if verbose:
			print "Not loading data (either it's been loaded already or it's not needed)"
	
	# And the simulation results
	if load_sim:
		if verbose:
			print "Loading simulation %s"%config["input"]["hoopoe"]
		hoopoe = s.shapecat(res=config["input"]["hoopoe"] ,truth=config["input"]["hoopoe"])
		hoopoe.res = fi.FITS(hoopoe.res_path)["i3s"].read()
		hoopoe.truth = fi.FITS(hoopoe.truth_path)["truth"].read()

		if (config["selection"]["weights"]):
			if (config["selection"]["weights_file"]!="none"):
				print "Reading weights"
				im3shape_weights = fi.FITS(config["selection"]["weights"])[-1].read()
				wts = im3shape_weights["weight"]
				hoopoe.res = arr.add_col(hoopoe.res, "weight", wts)
			else:
				print "Using im3shape weights from catalogue"
				print "mean weight : %3.3f"%(hoopoe.res["weight"].mean())

		sel = np.isfinite(hoopoe.res["mean_hsm_psf_e1_sky"]) & np.isfinite(hoopoe.res["mean_hsm_psf_e2_sky"])
		hoopoe.truth = hoopoe.truth[sel]
		hoopoe.res = hoopoe.res[sel]

		if (config["selection"]["mask"].lower()!="none"):
			apply_selection = True
			selection = fi.FITS(config["selection"]["mask"])["sel"].read().astype(bool)
			weights = fi.FITS(config["selection"]["mask"])["wts"].read()

			if verbose:
				print "Applying additional cuts and weights from %s"%config["selection"]["mask"]
			hoopoe.res = hoopoe.res[selection]
			hoopoe.truth = hoopoe.truth[selection]
			weights = weights[selection]
		if (not config["selection"]["reweight"]):
			if verbose:
				print "Ignoring weights."
			weights = np.ones(hoopoe.res["coadd_objects_id"].size)

		if not config["calibration"]["ztype"]:
			if verbose:
				print "Using DES redshift bins"

			exclude = (hoopoe.res["des_bin"]!=0 )
			hoopoe.truth = hoopoe.truth[exclude]  
			weights = weights[exclude]  
			hoopoe.res = hoopoe.res[exclude]
		else:
			if verbose:
				print "Using tophat redshift bins"

		if (config["selection"]["resample"]):
			print "Will apply resampling to match data"
			edat = np.sqrt(y1v2.res["e1"]**2+y1v2.res["e2"]**2)
			eh = np.sqrt(hoopoe.res["e1"]**2+hoopoe.res["e2"]**2)
			subsample = di.get_selection_to_match(edat,eh,nbins=35)
			hoopoe.res = hoopoe.res[subsample]
			hoopoe.truth = hoopoe.truth[subsample]
			weights = weights[subsample]


		print "Final selection : %d galaxies"%hoopoe.res["coadd_objects_id"].size
		print "Final selection : %d unique COSMOS IDs"%np.unique(hoopoe.truth["cosmos_ident"]).size

	else: 
		if verbose:
			print "Not loading simulation."

	return hoopoe, weights, y1v2


def choose_inputs(args):
	sim, data = False, False
	if args.calculate:
		sim = True
	if args.catalogue:
		data = True
	if (args.weights.lower()=="data"):
		data = True
	elif (args.weights.lower()=="sim"):
		sim = True

	return sim, data

def choose_outputs(config):
	"""Construct a list of outputs"""
	outputs = []
	for name in config["output"].keys():
		if (name=="dir") or (name=="filename"):
			continue
		else:
			if config["output"][name]:
				outputs.append(name)

	print "Will output the following :"
	print outputs
	return outputs

def main(args):

	load_sim, load_data = choose_inputs(args)
	outputs = choose_outputs(config)
	hoopoe, weights, y1v2 = setup(load_sim, load_data, config)
	
	if args.calculate:

		rbins= config["calibration"]["rbins"]
		sbins= config["calibration"]["sbins"]
		print "Using %d SNR bins , %d size bins"%(sbins,rbins)

		calculate(y1v2, hoopoe, 
			split_method=config["calibration"]["split"], 
			weights=weights,
			outputs=outputs,
			method=config["calibration"]["method"],
			config=config,
			sbins=sbins,
			rbins=rbins)
		

	if args.catalogue:
		rbins= config["calibration"]["rbins"]
		sbins= config["calibration"]["sbins"]
		print "Using %d SNR bins , %d size bins"%(sbins,rbins)
		calibrate(y1v2, config=config, sbins=sbins, rbins=rbins)

	if args.weights:
		hoopoe = s.shapecat()
		hoopoe.res = None
		get_weights(y1v2, hoopoe, config=config)

def mkdirs(config):
	os.system("mkdir -p %s/release/inputs"%config["output"]["dir"] )
	os.system("mkdir -p %s/release/outputs"%config["output"]["dir"] )
	os.system("mkdir -p %s/release/rbf"%config["output"]["dir"] )
	os.system("mkdir -p %s/release/grid"%config["output"]["dir"] )
	os.system("mkdir -p %s/release/polynomial"%config["output"]["dir"] )
	os.system("mkdir -p %s/nbc_data"%config["output"]["dir"] )


def calibrate(data, method="grid", config=None, smoothing=3, sbins=16, rbins=16):

	rbf = (method.lower()=="rbf")

	# Set up the wrappers
	nbc = cal.nbc()
	nbc_disc = cal.nbc()
	nbc.res = data.res
	del(data.res)
	gc.collect()
	nbc_bulge = cal.nbc()

	nbc_bulge.res=nbc.res[nbc.res["is_bulge"].astype(bool)]
	nbc_disc.res=nbc.res[np.invert(nbc.res["is_bulge"].astype(bool))]
	
	# Fit or interpolate
	if method is "polynomial":
		nbc_disc.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))
		nbc_disc.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))
		nbc_bulge.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))
		nbc_bulge.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))
	elif method is "rbf":
		nbc_disc.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins), smoothing=smoothing)
		nbc_bulge.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins), smoothing=smoothing)
	elif method is "grid":
		nbc_disc.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))[1].read()
		nbc_bulge.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))[1].read()

	# Apply to get a bias correction for each galaxy in the data
	nbc_disc.apply(split_half=0, scheme=method)
	print "Done disc calibration"

	nbc_bulge.apply(split_half=0, scheme=method)
	print "Done bulge calibration"

	print "Merging bulge/disc results...",
	nbc.res = arr.add_col(nbc.res, "m")
	nbc.res = arr.add_col(nbc.res, "c1")
	nbc.res = arr.add_col(nbc.res, "c2")
	nbc.res = arr.add_col(nbc.res, "a")

	nbc.res["m"][np.invert(nbc.res["is_bulge"].astype(bool))] = nbc_disc.res["m"]
	nbc.res["c1"][np.invert(nbc.res["is_bulge"].astype(bool))] = nbc_disc.res["c1"]
	nbc.res["c2"][np.invert(nbc.res["is_bulge"].astype(bool))] = nbc_disc.res["c2"]
	nbc.res["a"][np.invert(nbc.res["is_bulge"].astype(bool))] = nbc_disc.res["a"]
	nbc.res["m"][nbc.res["is_bulge"].astype(bool)] = nbc_bulge.res["m"]
	nbc.res["c1"][nbc.res["is_bulge"].astype(bool)] = nbc_bulge.res["c1"]
	nbc.res["c2"][nbc.res["is_bulge"].astype(bool)] = nbc_bulge.res["c2"]
	nbc.res["a"][nbc.res["is_bulge"].astype(bool)] = nbc_bulge.res["a"]
	print "done"


	if not os.path.exists(os.path.dirname(config["output"]["dir"])):
		os.system("mkdir -p %s"%os.path.dirname(config["output"]["dir"]))

	print "Saving calibrated catalogue to %s"%config["output_catalogue"]
	nbc.export(filename="%s/%s"%(config["output"]["dir"],config["output_catalogue"]))


def calculate(y1v2, hoopoe, split_method="none", weights=None, outputs=[], method="grid", smoothing=3, config=None, names=["m", "a"], sbins=16, rbins=16):

	if not os.path.exists(config["output"]["dir"]):
		print "Warning - output plots directory does not exist"

	# First make some diagnostic distributions of the observable parameters relevant for the NBC
	if ("histograms" in outputs):
		pl.histograms_vs_input(["e", "size", "flux"], hoopoe.truth, data2=y1v2.res, outdir="%s/release/inputs"%config["output"]["dir"], weights=weights)
		pl.histograms(["snr", "e", "size", "rgpp", "flux", "psfe", "psf_size"], hoopoe.res, kl=True, data2=y1v2.res, outdir="%s/release/outputs"%config["output"]["dir"], weights=weights)

	if ("alpha" in outputs):
		# Global fit of alpha to compare with previous runs
		plt.close()
		b = bias=di.get_alpha(hoopoe.res, hoopoe.res, nbins=15, names=["alpha11", "alpha22" ], xlim=(-0.03,0.02), xdata_name="mean_hsm_psf_e%d_sky", weights=weights, visual=True)
		plt.subplots_adjust(wspace=0, hspace=0, top=0.95, bottom=0.1)
		plt.title("Y1 v2 Sim, PSF v02")

		os.system("mkdir -p %s/release/alpha"%config["output"]["dir"])
		plt.savefig("%s/release/alpha/alphaplot-e-vs-epsf-linfit0.png"%config["output"]["dir"])
		plt.close()

	nbc = cal.nbc()
	nbc_disc = cal.nbc()
	nbc_bulge = cal.nbc()

	# split the data into two subsets
	wt1, wt2 = nbc.get_split_data(hoopoe, weights=weights, method=split_method)

	nbc_disc.load_from_cat(nbc, name="all")
	nbc_bulge.load_from_cat(nbc, name="all")

	#### Process disc then bulge runs

	# Fit the binned galaxies. Save one set of calibration data for disc objects, one for bulges
	edges_bulge = "equal_number"
	edges_disc = "equal_number"

	if ("half_tables" in outputs):
		nbc_disc.compute(split_half=1, 
			    fit="disc", weights=wt1,
			    use_catalogue_weights=config["selection"]["weights"], 
				reweight_per_bin=False, resample_per_bin=False,
				refdata=y1v2, binning=edges_disc, 
				rbins=rbins, sbins=sbins, 
				rlim=(1.13,3.0), slim=(12,200),
				table_name="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"], split_method, sbins,rbins))
		nbc_bulge.compute(split_half=1, 
			    fit="bulge", weights=wt1,
			    use_catalogue_weights=config["selection"]["weights"],
				reweight_per_bin=False, resample_per_bin=False,
				refdata=y1v2, binning=edges_bulge,
				rbins=rbins, sbins=sbins,
				rlim=(1.13,3.0), slim=(12,200),
				table_name="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"], split_method, sbins,rbins))

	if ("tables" in outputs):
		nbc_disc.compute(split_half=0,
			    fit="disc", weights=weights,
			    use_catalogue_weights=config["selection"]["weights"],
				reweight_per_bin=False, resample_per_bin=False,
				refdata=y1v2, binning=edges_disc,
				rbins=rbins, sbins=sbins,
				rlim=(1.13,3.0), slim=(12,200),
				table_name="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"], sbins, rbins))
		nbc_bulge.compute(split_half=0,
			    fit="bulge", weights=weights,
			    use_catalogue_weights=config["selection"]["weights"],
				reweight_per_bin=False, resample_per_bin=False,
				refdata=y1v2, binning=edges_bulge,
				rbins=rbins, sbins=sbins,
				rlim=(1.13,3.0), slim=(12,200),
				table_name="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"], sbins,rbins))


	sub_dir = method.lower()
	# Now plot out the points and the resulting smoothing fit
	if ("snr" in outputs):
		if (method.lower()=="rbf") or (method.lower()=="grid"):
			nbc_disc.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins), smoothing=smoothing)
			nbc_bulge.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins), smoothing=smoothing)
		elif (method.lower()=="polynomial"):
			nbc_disc.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))
			nbc_disc.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))
			nbc_bulge.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))
			nbc_bulge.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))

		nbc_disc.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins), bias_name="m", do_half=1, output="%s/m-vs-snr-disc-v1-1-s%2.2f-sbins%d-rbins%d.png"%(config["output"]["dir"]+"/release/"+sub_dir, smoothing, sbins, rbins), use_rbf=rbf)
		nbc_disc.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins), bias_name="a", do_half=1, output="%s/alpha-vs-snr-disc-v1-1.png"%(config["output"]["dir"]+"/release/"+ sub_dir), use_rbf=rbf)
		plt.close()
		nbc_disc.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins), bias_name="m", do_half=2, output="%s/m-vs-snr-disc-v1-2-s%2.2f-sbins%d-rbins%d.png"%(config["output"]["dir"]+"/release/"+ sub_dir, smoothing, sbins, rbins), use_rbf=rbf)
		nbc_disc.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins), bias_name="a", do_half=2, output="%s/alpha-vs-snr-disc-v1-2.png"%(config["output"]["dir"]+"/release/"+ sub_dir), use_rbf=rbf)
		plt.close()
		nbc_bulge.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins), bias_name="m", do_half=1, output="%s/m-vs-snr-bulge-v1-2-s%2.2f-sbins%d-rbins%d.png"%(config["output"]["dir"]+"/release/"+ sub_dir, smoothing, sbins, rbins), use_rbf=rbf)
		nbc_bulge.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins), bias_name="a", do_half=1, output="%s/alpha-vs-snr-bulge-v1-1.png"%(config["output"]["dir"]+"/release/"+ sub_dir), use_rbf=rbf)
		plt.close()
		nbc_bulge.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins), bias_name="m", do_half=2, output="%s/m-vs-snr-bulge-v1-1-s%2.2f-sbins%d-rbins%d.png"%(config["output"]["dir"]+"/release/"+ sub_dir, smoothing, sbins, rbins), use_rbf=rbf)
		nbc_bulge.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins), bias_name="a", do_half=2, output="%s/alpha-vs-snr-bulge-v1-2.png"%(config["output"]["dir"]+"/release/"+ sub_dir), use_rbf=rbf)
		plt.close()

	if (split_method!="none"):
		# Derive a calibration from the grid in whichever way is specified in the config file
		if (method=="rbf"):
			nbc_disc.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"], split_method,sbins,rbins), smoothing=smoothing)
			nbc_bulge.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"], split_method,sbins,rbins), smoothing=smoothing)
		elif (method=="grid"):
			nbc_disc.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"], split_method,sbins,rbins))[1].read()
			nbc_bulge.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"], split_method,sbins,rbins))[1].read()
		elif (method=="polynomial"):
			nbc_disc.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"], split_method,sbins,rbins))
			nbc_disc.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"], split_method,sbins,rbins))
			nbc_bulge.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"], split_method,sbins,rbins))
			nbc_bulge.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"], split_method,sbins,rbins))

		# Apply it to the other half
		nbc_disc.apply(split_half=2, scheme=method, names=names)
		nbc_bulge.apply(split_half=2, scheme=method, names=names)
		nbc.combine_bd(nbc_disc,nbc_bulge, split_half=2, names=["m"]*("m" in names) + ["c1", "c2"]*("a" in names) )

	# Finally save some diagnostic plots in tomographic bins
	if ("redshift" in outputs): 
		zbins=[ 0.2, 0.43, 0.63, 0.9, 1.3]
		tophat = config["calibration"]["ztype"]=="tophat"

		import pdb ; pdb.set_trace()

		# Only do the redshift plots with the split catalogues if they were recalculated this time
		if ("half_tables" in outputs):
			plt.close()
			if "m" in names:
				bias0=nbc.redshift_diagnostic(bias="m", label="Uncalibrated", ls="none", nbins=4, fmt=["o","D"], colour="steelblue", weights=wt2, bins=zbins, tophat=tophat, separate_components=False)
				bias=nbc.redshift_diagnostic(bias="m", label="Calibrated", ls="none", nbins=4, fmt=["^",">"], apply_calibration=True, colour="purple", weights=wt2, bins=zbins, tophat=tophat, separate_components=False)
				plt.ylabel("Multiplicative Bias $m$")
				plt.legend(loc="center right")
				plt.savefig("%s/release/%s/m-bias-vs-redshift-diagnostic-v1-%s-halfcat-s%2.3f-sbins%d-rbins%d-tophat%d.png"%(config["output"]["dir"], sub_dir, split_method, smoothing, sbins,rbins, int(tophat)))
				plt.close()
			if "a" in names:
				nbc.redshift_diagnostic(bias="alpha", label="Uncalibrated",ls="none", nbins=4, fmt=["o","D"], colour="steelblue", weights=wt2, bins=zbins, tophat=tophat)
				nbc.redshift_diagnostic(bias="alpha", label="Calibrated", ls="none",fmt=["^",">"], nbins=4, apply_calibration=True, colour="purple", weights=wt2, bins=zbins, tophat=tophat)
				plt.ylabel(r"PSF Leakage $\alpha$")
				plt.legend(loc="upper left")
				plt.savefig("%s/release/%s/alpha-vs-redshift-diagnostic-v1-%s-halfcat-s%2.3f-sbins%d-rbins%d-tophat%d.png"%(config["output"]["dir"], sub_dir, split_method, smoothing, sbins,rbins, int(tophat)))
				plt.close()

		# Finally redo the fits with the full catalogue
		nbc = cal.nbc()
		nbc_disc = cal.nbc()
		nbc_bulge = cal.nbc()
		wt1, wt2 = nbc.get_split_data(hoopoe, weights=weights)

		nbc_disc.load_from_cat(nbc, name="all")
		nbc_bulge.load_from_cat(nbc, name="all")

		if (method.lower()=="rbf"):
			nbc_disc.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins), smoothing=smoothing)
			nbc_bulge.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins), smoothing=smoothing)
		elif (method.lower()=="grid"):
			nbc_disc.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))[1].read()
			nbc_bulge.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))[1].read()
		elif (method.lower()=="polynomial"):
			nbc_disc.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))
			nbc_disc.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))
			nbc_bulge.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))
			nbc_bulge.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output"]["dir"],sbins,rbins))


		nbc_disc.apply(split_half=0, scheme=method, names=names)
		nbc_bulge.apply(split_half=0, scheme=method, names=names)
		nbc.combine_bd(nbc_disc,nbc_bulge, split_half=0, names=["m"]*("m" in names) + ["c1", "c2"]*("a" in names) )

		plt.close()
		if "m" in names:
			bias0 = nbc.redshift_diagnostic(bias="m", label="Uncalibrated", ls="none", nbins=4, fmt=["o","D"], colour="steelblue", weights=weights, split_half=0, bins=zbins, tophat=tophat, separate_components=False)
			bias = nbc.redshift_diagnostic(bias="m", label="Calibrated", ls="none", nbins=4, fmt=["^",">"], apply_calibration=True, colour="purple", weights=weights, split_half=0, bins=zbins, tophat=tophat, separate_components=False)
			plt.ylabel("Multiplicative Bias $m$")
			plt.legend(loc="center right")
			plt.savefig("%s/release/%s/m-bias-vs-redshift-diagnostic-v1-fullcat-s%2.2f-sbins%d-rbins%d-tophat%d.png"%(config["output"]["dir"],sub_dir, smoothing, sbins, rbins, int(tophat)))
			plt.close()

		if "a" in names:
			nbc.redshift_diagnostic(bias="alpha", label="Uncalibrated",ls="none", nbins=4, fmt=["o","D"], colour="steelblue", weights=weights, split_half=0, bins=zbins, tophat=tophat)
			nbc.redshift_diagnostic(bias="alpha", label="Calibrated", ls="none",fmt=["^",">"], nbins=4, apply_calibration=True, colour="purple", weights=weights, split_half=0, bins=zbins, tophat=tophat)
			plt.ylabel(r"PSF Leakage $\alpha$")
			plt.legend(loc="upper left")
			plt.savefig("%s/release/%s/alpha-vs-redshift-diagnostic-v1-fullcat-s%2.2f-sbins%d-rbins%d-tophat%d.png"%(config["output"]["dir"],sub_dir, smoothing, sbins, rbins, int(tophat)))
			plt.close()


def get_weights(y1v2, hoopoe, sbins=9, rbins=9, config=None, catalogue="data"):
	os.system("mkdir -p %s/weights/"%config["output"]["dir"])
	import pdb ; pdb.set_trace()
	weights_grid = di.im3shape_weights_grid(y1v2.res, bins_from_table=False, filename="%s/weights/%sim3shape_weights_grid_v5_extra.fits"%(config["output"]["dir"],catalogue), sbins=sbins, rbins=rbins, simdat=hoopoe.res, binning="log")
	di.interpolate_weights_grid(weights_grid, y1v2.res, smoothing=2.5, outdir="%s/weights/"%config["output"]["dir"], outfile="%s-hoopoe_weights_column-v4_extra.fits"%catalogue)



def bin_edges_from_table(table_dir, type="disc"):
	table_name = glob.glob("%s/bias_table_hoopoe-v1-fullcat-%s.fits"%(table_dir,type))
	if len(table_name)==0:
		print "No table of name %s/bias_table_hoopoe-v1-fullcat-%s.fits"%(table_dir,type)
	else:
		table_name = table_name[0]
		print "Reading fixed bin edges from %s"%table_name

	table = fi.FITS(table_name)[1].read()

	rgpp_edges = zip(np.log10(table["rgp_lower"]), np.log10(table["rgp_upper"]))
	snr_edges = zip(np.log10(table["snr_lower"]), np.log10(table["snr_upper"]))

	print rgpp_edges
	print snr_edges

	return rgpp_edges, snr_edges



if __name__ == "__main__":
	parser = argparse.ArgumentParser(add_help=False)
	parser.add_argument('--config',"-c", type=str, action='store')
	parser.add_argument('--catalogue', action='store_true')
	parser.add_argument('--calculate', action='store_true')
	parser.add_argument('--weights', type=str, action='store', default='none')

	args = parser.parse_args()

	config =yaml.load(open(args.config))

	mkdirs(config)
	main(args)

