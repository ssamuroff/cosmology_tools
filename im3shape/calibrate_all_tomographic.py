import numpy as np
import tools.nbc as cal
import tools.shapes as s
import tools.diagnostics as di
import tools.plots as pl
import tools.arrays as arr
import fitsio as fi
import pylab as plt
import os, yaml, argparse, glob, gc, copy

plt.switch_backend("agg")


def setup(calculate, catalogue, config, y1v2=None):
	im3shape_columns = ["e1", "e2", "mean_hsm_psf_e1_sky", "mean_hsm_psf_e2_sky", "mean_hsm_psf_sigma", "snr", "mean_rgpp_rp", "coadd_objects_id", "is_bulge", "bulge_flux", "disc_flux", "mean_psf_e1_sky", "mean_psf_e2_sky", "mean_psf_fwhm", "info_flag"]
	truth_columns = ['DES_id', 'cosmos_ident', 'cosmos_photoz', 'sextractor_pixel_offset', 'true_g1', 'true_g2', 'intrinsic_e1', 'intrinsic_e2', 'ra', 'dec', 'hlr', 'mag', 'flux']

	# Load the y1 data
	if (y1v2 is  None ) and ((calculate and config["output"]["histograms"]) or catalogue or config["resample"] or config["resample_perbin"] or config["reweight_perbin"]):
		y1v2 = s.shapecat(res=config["i3s_dir"])
		y1v2.load(truth=False, prune=True, cols=[im3shape_columns,truth_columns])
		#y1v2.res=y1v2.res[y1v2.res["info_flag"]==0]
		sel=((y1v2.res["snr"] > 12) & (y1v2.res["snr"] < 200) & (y1v2.res["mean_rgpp_rp"] > 1.13) & (y1v2.res["mean_rgpp_rp"] < 3.0))
		y1v2.res=y1v2.res[sel]
	else:
		print "Not loading data (either it's been loaded already or it's not needed)"
	

	# And the simulation results
	if calculate:
		hoopoe = s.shapecat(res=config["hoopoe_dir"] ,truth=config["hoopoe_dir"])
		hoopoe.res=fi.FITS(hoopoe.res_path)["i3s"].read()
		hoopoe.truth=fi.FITS(hoopoe.truth_path)["truth"].read()

		sel = np.isfinite(hoopoe.res["mean_hsm_psf_e1_sky"]) & np.isfinite(hoopoe.res["mean_hsm_psf_e2_sky"])
		hoopoe.truth = hoopoe.truth[sel]
		hoopoe.res = hoopoe.res[sel]
		

		apply_selection = False
		if ("apply_selection" in config.keys()):
			if config["apply_selection"]:
				apply_selection = True
				selection = fi.FITS(config["selection"])["sel"].read().astype(bool)
				weights = fi.FITS(config["selection"])["wts"].read()

				print "Applying additional cuts and weights from %s"%config["selection"]
				hoopoe.res = hoopoe.res[selection]
				hoopoe.truth = hoopoe.truth[selection]
				weights = weights[selection]

		if (not apply_selection) or (not config["reweight"]):
			weights = np.ones(hoopoe.res["coadd_objects_id"].size)

		if not config["tophat_binning"]:
			print "Using DES redshift bins"

			#bin_allocation = fi.FITS("/share/des/disc6/samuroff/y1/hoopoe/hoopoe-v2-zbin_allocation.fits")[1].read()
			#bin_num, hoopoe.res = di.match_results(bin_allocation, hoopoe.res)
			#hoopoe.res = arr.add_col(hoopoe.res,"des_bin", bin_num["bin"])
			exclude = (hoopoe.res["des_bin"]!=0 )
			hoopoe.truth = hoopoe.truth[exclude]  
			weights = weights[exclude]  
			hoopoe.res = hoopoe.res[exclude]
		else:
			print "Using tophat redshift bins"

		if (config["resample"]):
			print "Will apply resampling to match data"
			edat = np.sqrt(y1v2.res["e1"]**2+y1v2.res["e2"]**2)
			eh = np.sqrt(hoopoe.res["e1"]**2+hoopoe.res["e2"]**2)
			subsample = di.get_selection_to_match(edat,eh,nbins=35)
			hoopoe.res = hoopoe.res[subsample]
			hoopoe.truth = hoopoe.truth[subsample]
			weights = weights[subsample]

		print "Final selection : %d galaxies"%hoopoe.res["coadd_objects_id"].size
		print "Final selection : %d unique COSMOS IDs"%np.unique(hoopoe.truth["cosmos_ident"]).size

		return hoopoe, weights, y1v2


def main(args):
	hoopoe, weights, y1v2 = setup(args.calculate, args.catalogue, config)
	
	# And the simulation results
	if args.calculate:
		rbins=11
		sbins=11

		diagnostics(y1v2, hoopoe, weights=weights, tomographic_calibration=config["tomography"], histograms=False, alpha=False, table=config["output"]["tables"], vsredshift=True, rbf=False, simple_grid=True, config=config, sbins=sbins, rbins=rbins, half_tables=config["random_halves"])

		diagnostics(y1v2, hoopoe, weights=weights, tomographic_calibration=config["tomography"], vssnr=config["output"]["snr"], vsredshift=config["output"]["redshift"], table=config["output"]["tables"], alpha=config["output"]["alpha"], histograms=config["output"]["histograms"], rbf=True, config=config, sbins=sbins, rbins=rbins, half_tables=config["random_halves"])
		diagnostics(y1v2, hoopoe, weights=weights, tomographic_calibration=config["tomography"], histograms=False, alpha=False, table=False, vsredshift=True, rbf=False, simple_grid=False, config=config, sbins=sbins, rbins=rbins, half_tables=config["random_halves"])
		

		if config["cosmos_halves"]:
			diagnostics(y1v2, hoopoe, split_method="cosmos", weights=weights, tomographic_calibration=config["tomography"], histograms=False, alpha=False, table=False, half_tables=config["output"]["tables"], vsredshift=config["output"]["redshift"], rbf=True, config=config, sbins=sbins, rbins=rbins)
			diagnostics(y1v2, hoopoe, split_method="cosmos", weights=weights, tomographic_calibration=config["tomography"], histograms=False, alpha=False, table=False, half_tables=False, vsredshift=config["output"]["redshift"], rbf=False, simple_grid=False, config=config, sbins=sbins, rbins=rbins)
			diagnostics(y1v2, hoopoe, split_method="cosmos", weights=weights, tomographic_calibration=config["tomography"], histograms=False, alpha=False, table=False, half_tables=False, vsredshift=config["output"]["redshift"], rbf=False, simple_grid=True, config=config, sbins=sbins, rbins=rbins)

	if args.catalogue:
		calibrate(y1v2)

	if args.weights:
		get_weights(y1v2, hoopoe, config=config)

def mkdirs(config):
	os.system("mkdir -p %s/release/inputs"%config["output_dir"] )
	os.system("mkdir -p %s/release/outputs"%config["output_dir"] )
	os.system("mkdir -p %s/release/rbf"%config["output_dir"] )
	os.system("mkdir -p %s/release/grid"%config["output_dir"] )
	os.system("mkdir -p %s/release/polynomial"%config["output_dir"] )
	os.system("mkdir -p %s/nbc_data"%config["output_dir"] )


def calibrate(data, method="rbf", smoothing=3, sbins=16, rbins=10):

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
		nbc_disc.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output_dir"],sbins,rbins))
		nbc_disc.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output_dir"],sbins,rbins))
		nbc_bulge.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output_dir"],sbins,rbins))
		nbc_bulge.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output_dir"],sbins,rbins))
	elif method is "rbf":
		nbc_disc.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output_dir"],sbins,rbins), smoothing=smoothing)
		nbc_bulge.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output_dir"],sbins,rbins), smoothing=smoothing)
	elif method is "grid":
		nbc_disc.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins.fits"%(config["output_dir"],sbins,rbins))[1].read()
		nbc_bulge.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins.fits"%(config["output_dir"],sbins,rbins))[1].read()

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


	if not os.path.exists(os.path.dirname(config["output_dir"])):
		os.system("mkdir -p %s"%os.path.dirname(config["output_dir"]))

	print "Saving calibrated catalogue to %s"%config["output_catalogue"]
	nbc.export(filename="%s/%s"%(config["output_dir"],config["output_catalogue"]))


def diagnostics(y1v2, hoopoe, tomographic_calibration=False, histograms=True, split_method="random", weights=None, alpha=True, table=True, half_tables=True, vssnr=True, vsredshift=True, rbf=True, simple_grid=False, smoothing=3, config=None, names=["m", "a"], sbins=16, rbins=11):

	if rbf:
		sub_dir="rbf"
	elif simple_grid:
		sub_dir="grid"
	else:
		sub_dir="polynomial"

	if not os.path.exists(config["output_dir"]):
		print "Warning - output plots directory does not exist"

	# First make some diagnostic distributions of the observable parameters relevant for the NBC
	if histograms:
		pl.histograms_vs_input(["e", "size", "flux"], hoopoe.truth, data2=y1v2.res, outdir="%s/release/inputs"%config["output_dir"], weights=weights)
		pl.histograms(["snr", "e", "size", "rgpp", "flux", "psfe", "psf_size"], hoopoe.res, kl=True, data2=y1v2.res, outdir="%s/release/outputs"%config["output_dir"], weights=weights)

	if alpha:
		# Global fit of alpha to compare with previous runs
		plt.close()
		b = bias=di.get_alpha(hoopoe.res, hoopoe.res, nbins=15, names=["alpha11", "alpha22" ], xlim=(-0.03,0.02), xdata_name="mean_hsm_psf_e%d_sky", weights=weights, visual=True)
		plt.subplots_adjust(wspace=0, hspace=0, top=0.95, bottom=0.1)
		plt.title("Y1 v2 Sim, PSF v02")

		os.system("mkdir -p %s/release/alpha"%config["output_dir"])
		plt.savefig("%s/release/alpha/alphaplot-e-vs-epsf-linfit0.png"%config["output_dir"])
		plt.close()


	#Exit here if no calibration diagnostics are needed
	if (not table) and (not vssnr) and (not vsredshift):
		return 0

	# Now compute the calibration using a random half of the simulation
	if vssnr or half_tables: 
		nbc = cal.nbc()
		nbc_disc = cal.nbc()
		nbc_bulge = cal.nbc()
		wt1, wt2 = nbc.get_split_data(hoopoe, weights=weights, method=split_method)

		nbc_disc.load_from_cat(nbc, name="all")
		nbc_bulge.load_from_cat(nbc, name="all")

	#### Process disc then bulge runs

	# Fit the binned galaxies. Save one set of calibration data for disc objects, one for bulges
	import pdb ; pdb.set_trace()
	if table:
		if ("match_gridpoints" in config.keys()):
			redges_bulge, sedges_bulge = bin_edges_from_table(config["match_gridpoints"], type="bulge")
			redges_disc, sedges_disc = bin_edges_from_table(config["match_gridpoints"], type="disc")

			edges_disc = (redges_disc, sedges_disc)
			edges_bulge = (redges_bulge, sedges_bulge)

		else:
			edges_bulge = "equal_number"
			edges_disc = "equal_number"

		if tomographic_calibration:
			zbin_numbers = np.unique(nbc.res["des_bin"])

		if half_tables:
			print "Calculating half tables in %d redshift bins"%zbin_numbers.size
			for k in zbin_numbers:
				print "bin %d"%k
				nbc_disc.compute(split_half=1, fit="disc", weights=wt1, reweight_per_bin=config["reweight_perbin"], resample_per_bin=config["resample_perbin"], refdata=y1v2, binning=edges_disc, redshift_bin=k, rbins=rbins, sbins=sbins, rlim=(1.13,3.0), slim=(12,200), table_name="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-disc-%dsbins-%drbins-zbin%d.fits"%(config["output_dir"], split_method, sbins,rbins, k))
				nbc_bulge.compute(split_half=1, fit="bulge", weights=wt1, reweight_per_bin=config["reweight_perbin"], resample_per_bin=config["resample_perbin"], refdata=y1v2, binning=edges_bulge, redshift_bin=k, rbins=rbins, sbins=sbins, rlim=(1.13,3.0), slim=(12,200), table_name="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-bulge-%dsbins-%drbins-zbin%d.fits"%(config["output_dir"], split_method, sbins,rbins, k))

		if table:
			print "Calculating full tables in %d redshift bins"%zbin_numbers.size
			for k in zbin_numbers:
				print "bin %d"%k
				nbc_disc.compute(split_half=0, fit="disc", weights=weights, reweight_per_bin=config["reweight_perbin"], resample_per_bin=config["resample_perbin"], refdata=y1v2, binning=edges_disc, redshift_bin=k, rbins=rbins, sbins=sbins, rlim=(1.13,3.0), slim=(12,200), table_name="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc-%dsbins-%drbins-zbin%d.fits"%(config["output_dir"], sbins, rbins, k))
				nbc_bulge.compute(split_half=0, fit="bulge", weights=weights, reweight_per_bin=config["reweight_perbin"], resample_per_bin=config["resample_perbin"], refdata=y1v2, binning=edges_bulge, redshift_bin=k, rbins=rbins, sbins=sbins, rlim=(1.13,3.0), slim=(12,200), table_name="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge-%dsbins-%drbins-zbin%d.fits"%(config["output_dir"], sbins,rbins, k))

	if rbf:
		scheme="rbf"
	elif simple_grid:
		scheme="grid"
	else:
		scheme="polynomial"

	
	if half_tables:
		nbc = calibrate_in_bins(nbc, halfcats=True, split_method=split_method, scheme=scheme, config=config, sbins=sbins, rbins=rbins, smoothing=smoothing, names=names)
			
	# Finally save some diagnostic plots in tomographic bins
	if vsredshift: 
		zbins=[ 0.2, 0.43, 0.63, 0.9, 1.3]
		tophat = config["tophat_binning"]

		if half_tables:
			plt.close()
			if "m" in names:
				nbc.redshift_diagnostic(bias="m", label="Uncalibrated", ls="none", nbins=3, fmt=["o","D"], colour="steelblue", weights=wt2, bins=zbins, tophat=tophat, separate_components=False)
				nbc.redshift_diagnostic(bias="m", label="Calibrated", ls="none", nbins=3, fmt=["^",">"], apply_calibration=True, colour="purple", weights=wt2, bins=zbins, tophat=tophat, separate_components=False)
				plt.ylabel("Multiplicative Bias $m$")
				plt.legend(loc="center right")
				plt.savefig("%s/release/%s/m-bias-vs-redshift-diagnostic-v1-%s-halfcat-s%2.3f-sbins%d-rbins%d-tophat%d-tomographic_calibration.png"%(config["output_dir"], sub_dir, split_method, smoothing, sbins,rbins, int(tophat)))
				plt.close()
			if "a" in names:
				nbc.redshift_diagnostic(bias="alpha", label="Uncalibrated",ls="none", nbins=3, fmt=["o","D"], colour="steelblue", weights=wt2, bins=zbins, tophat=tophat)
				nbc.redshift_diagnostic(bias="alpha", label="Calibrated", ls="none",fmt=["^",">"], nbins=3, apply_calibration=True, colour="purple", weights=wt2, bins=zbins, tophat=tophat)
				plt.ylabel(r"PSF Leakage $\alpha$")
				plt.legend(loc="upper left")
				plt.savefig("%s/release/%s/alpha-vs-redshift-diagnostic-v1-%s-halfcat-s%2.3f-sbins%d-rbins%d-tophat%d-tomographic_calibration.png"%(config["output_dir"], sub_dir, split_method, smoothing, sbins,rbins, int(tophat)))
				plt.close()

		# Finally redo the fits with the full catalogue
		gc.collect()
		nbc = cal.nbc()
		wt1, wt2 = nbc.get_split_data(hoopoe, weights=weights, method=split_method)
		nbc = calibrate_in_bins(nbc, halfcats=False, split_method=None, scheme=scheme, config=config, sbins=sbins, rbins=rbins, smoothing=smoothing, names=names)

		import pdb ; pdb.set_trace()

		plt.close()
		if "m" in names:
			nbc.redshift_diagnostic(bias="m", label="Uncalibrated", ls="none", nbins=3, fmt=["o","D"], colour="steelblue", weights=weights, split_half=0, bins=zbins, tophat=tophat, separate_components=False)
			bias=nbc.redshift_diagnostic(bias="m", label="Calibrated", ls="none", nbins=3, fmt=["^",">"], apply_calibration=True, colour="purple", weights=weights, split_half=0, bins=zbins, tophat=tophat, separate_components=False)
			plt.ylabel("Multiplicative Bias $m$")
			plt.legend(loc="center right")
			plt.savefig("%s/release/%s/m-bias-vs-redshift-diagnostic-v1-fullcat-s%2.2f-sbins%d-rbins%d-tophat%d-tomographic_calibration.png"%(config["output_dir"],sub_dir, smoothing, sbins, rbins, int(tophat)))
			plt.close()

		if "a" in names:
			nbc.redshift_diagnostic(bias="alpha", label="Uncalibrated",ls="none", nbins=3, fmt=["o","D"], colour="steelblue", weights=weights, split_half=0, bins=zbins, tophat=tophat)
			nbc.redshift_diagnostic(bias="alpha", label="Calibrated", ls="none",fmt=["^",">"], nbins=3, apply_calibration=True, colour="purple", weights=weights, split_half=0, bins=zbins, tophat=tophat)
			plt.ylabel(r"PSF Leakage $\alpha$")
			plt.legend(loc="upper left")
			plt.savefig("%s/release/%s/alpha-vs-redshift-diagnostic-v1-fullcat-s%2.2f-sbins%d-rbins%d-tophat%d-tomographic_calibration.png"%(config["output_dir"],sub_dir, smoothing, sbins, rbins, int(tophat)))
			plt.close()


def get_weights(y1v2, hoopoe, sbins=21, rbins=21, config=None):
	os.system("mkdir -p %s/weights/"%config["output_dir"])
	weights_grid = di.im3shape_weights_grid(y1v2, bins_from_table=False, filename="%s/weights/im3shape_weights_grid.fits"%config["output_dir"], sbins=sbins, rbins=rbins, simdat=hoopoe.res, binning="fixed")
	di.interpolate_weights_grid(weights_grid, y1v2, smoothing=1, outdir="%s/weights/"%config["output_dir"], outfile="im3shape_weights_column.fits")



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

def set_defaults(config):
	names=["resample", "resample_perbin", "reweight", "reweight_perbin", "tophat_binning"]
	defaults = [False,False,False,False,False]
	for (name,default) in zip(names,defaults):
		if name not in config.keys(): config[name] = default

	return config


def calibrate_in_bins(nbc, halfcats=False, split_method="", scheme="rbf", config=None, sbins=12, rbins=9, smoothing=3, names=["m", "a"]):

	discs = cal.nbc()
	bulges = cal.nbc()

	if halfcats: 
		split_half = 2
		prefix = "%s-halfcat"%split_method
	else:
		split_half = 0
		prefix = "fullcat"

	zbin_numbers = np.unique(nbc.res["des_bin"])
	print "Calibrating in %d redshift bins"%zbin_numbers.size
	for k in zbin_numbers:
		print "-- calibrating galaxies in bin %d"%k
		

		if halfcats:
			select_bin1 = nbc.res1["des_bin"]==k
			select_bin = nbc.res2["des_bin"]==k
			discs.res1 = copy.deepcopy(nbc.res1[select_bin1])
			discs.truth1 = copy.deepcopy(nbc.truth1[select_bin1])
			bulges.res1 = copy.deepcopy(nbc.res1[select_bin1])
			bulges.truth1 = copy.deepcopy(nbc.truth1[select_bin1])
			discs.res2 = copy.deepcopy(nbc.res2[select_bin])
			discs.truth2 = copy.deepcopy(nbc.truth2[select_bin])
			bulges.res2 = copy.deepcopy(nbc.res2[select_bin])
			bulges.truth2 = copy.deepcopy(nbc.truth2[select_bin])

		else:
			select_bin = nbc.res["des_bin"]==k
			discs.res = copy.deepcopy(nbc.res[select_bin])
			discs.truth = copy.deepcopy(nbc.truth[select_bin])
			bulges.res = copy.deepcopy(nbc.res[select_bin])
			bulges.truth = copy.deepcopy(nbc.truth[select_bin])

		if (scheme=="rbf"):
			discs.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-%s-disc-%dsbins-%drbins-zbin%d.fits"%(config["output_dir"], prefix,sbins,rbins, k), smoothing=smoothing)
			bulges.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-%s-bulge-%dsbins-%drbins-zbin%d.fits"%(config["output_dir"], prefix,sbins,rbins, k), smoothing=smoothing)
		elif (scheme=="grid"):
			discs.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-%s-disc-%dsbins-%drbins-zbin%d.fits"%(config["output_dir"], prefix,sbins,rbins, k))[1].read()
			bulges.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-%s-bulge-%dsbins-%drbins-zbin%d.fits"%(config["output_dir"], prefix,sbins,rbins, k))[1].read()
		elif (scheme=="polynomial"):
			discs.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-%s-disc-%dsbins-%drbins-zbin%d.fits"%(config["output_dir"], prefix,sbins,rbins, k))
			discs.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-%s-disc-%dsbins-%drbins-zbin%d.fits"%(config["output_dir"], prefix,sbins,rbins, k))
			bulges.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-%s-bulge-%dsbins-%drbins-zbin%d.fits"%(config["output_dir"], prefix,sbins,rbins, k))
			bulges.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-%s-bulge-%dsbins-%drbins-zbin%d.fits"%(config["output_dir"], prefix,sbins,rbins, k))

		discs.apply(split_half=split_half, scheme=scheme, names=names)
		bulges.apply(split_half=split_half, scheme=scheme, names=names)
		nbc.combine_bd(discs,bulges, split_half=split_half, mask=select_bin, names=["m"]*("m" in names) + ["c1", "c2"]*("a" in names) )

	return nbc

def generate_difference_table(config1, config2,sbins=16, rbins=16, y1v2=None):
	hoopoe1, weights1, y1v2 = ct.setup(True, True, config1, y1v2=y1v2)
	hoopoe2, weights2, y1v2 = ct.setup(True, False, config2, y1v2=y1v2)

	subset = np.in1d(hoopoe1.res["coadd_objects_id"], hoopoe2.res["coadd_objects_id"])

	nbc = cal.nbc()
	nbc_disc = cal.nbc()
	nbc_bulge = cal.nbc()
	wt1, wt2 = nbc.get_split_data(hoopoe1, weights=weights1, method="random")
	nbc_disc.load_from_cat(nbc, name="all")
	nbc_bulge.load_from_cat(nbc, name="all")

	import pdb ; pdb.set_trace()

	nbc_disc.diff_table(weights=weights1, compare_with_subset=subset, fit="disc", table_name="diff_table_disc.fits", sbins=sbins, rbins=rbins, binning="equal_number", rlim=(1.13,3), slim=(12,200))
	nbc_bulge.diff_table(weights=weights2, compare_with_subset=subset, fit="bulge", table_name="diff_table_bulge.fits", sbins=sbins, rbins=rbins, binning="equal_number", rlim=(1.13,3), slim=(12,200))

if __name__ == "__main__":
	parser = argparse.ArgumentParser(add_help=False)
	parser.add_argument('--config',"-c", type=str, action='store')
	parser.add_argument('--catalogue', action='store_true')
	parser.add_argument('--calculate', action='store_true')
	parser.add_argument('--weights', action='store_true')

	args = parser.parse_args()

	config = set_defaults(yaml.load(open(args.config)))

	mkdirs(config)
	main(args)

