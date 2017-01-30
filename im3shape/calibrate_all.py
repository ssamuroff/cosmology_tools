import numpy as np
import tools.nbc as cal
import tools.shapes as s
import tools.diagnostics as di
import tools.plots as pl
import fitsio as fi
import pylab as plt
import os, yaml, argparse

plt.switch_backend("agg")


def main(args):
	im3shape_columns = ["e1", "e2", "mean_hsm_psf_e1_sky", "mean_hsm_psf_e2_sky", "mean_hsm_psf_sigma", "snr", "mean_rgpp_rp", "radius", "coadd_objects_id", "mean_flux", "n_exposure", "stamp_size", "info_flag", "is_bulge", "tilename"]
	truth_columns = ['DES_id', 'cosmos_ident', 'cosmos_photoz', 'sextractor_pixel_offset', 'true_g1', 'true_g2', 'intrinsic_e1', 'intrinsic_e2', 'ra', 'dec', 'hlr', 'mag', 'flux']

	# Load the y1 data
	if args.calculate or args.catalogue:
		y1v2 = s.shapecat(res=config["i3s_dir"])
		y1v2.load(truth=False, prune=True, cols=[im3shape_columns,truth_columns])
		y1v2.res=y1v2.res[y1v2.res["info_flag"]==0]
		sel=((y1v2.res["snr"] > 12) & (y1v2.res["snr"] < 200) & (y1v2.res["mean_rgpp_rp"] > 1.13) & (y1v2.res["mean_rgpp_rp"] < 3.0))
		y1v2.res=y1v2.res[sel]
	

	# And the simulation results
	if args.calculate:
		hoopoe = s.shapecat(res=config["hoopoe_dir"] ,truth=config["hoopoe_dir"])
		hoopoe.res=fi.FITS(hoopoe.res_path)["i3s"].read()
		hoopoe.truth=fi.FITS(hoopoe.truth_path)["truth"].read()
		sel = np.isfinite(hoopoe.res["mean_hsm_psf_e1_sky"]) & np.isfinite(hoopoe.res["mean_hsm_psf_e2_sky"])
		hoopoe.res = hoopoe.res[sel]
		hoopoe.truth = hoopoe.truth[sel]

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

				print "Final selection : %d galaxies"%hoopoe.res["coadd_objects_id"].size
				print "Final selection : %d unique COSMOS IDs"%np.unique(hoopoe.truth["cosmos_ident"]).size

		if (not apply_selection) or (not config["reweight"]):
			weights = np.ones(hoopoe.res["coadd_objects_id"].size)

		if (config["resample"]):
			edat = np.sqrt(y1v2.res["e1"]**2+y1v2.res["e2"]**2)
			eh = np.sqrt(hoopoe.res["e1"]**2+hoopoe.res["e2"]**2)
			subsample = di.get_selection_to_match(edat,eh,nbins=35)
			hoopoe.res = hoopoe.res[subsample]
			hoopoe.truth = hoopoe.truth[subsample]

		print "Final selection : %d galaxies"%hoopoe.res["coadd_objects_id"].size
		print "Final selection : %d unique COSMOS IDs"%np.unique(hoopoe.truth["cosmos_ident"]).size

		diagnostics(y1v2, hoopoe, weights=weights, vssnr=config["output"]["snr"], vsredshift=config["output"]["redshift"], table=config["output"]["tables"], alpha=config["output"]["alpha"], histograms=config["output"]["histograms"], rbf=True)
		diagnostics(y1v2, hoopoe, weights=weights, histograms=False, alpha=False, table=False, rbf=False, simple_grid=False)
		diagnostics(y1v2, hoopoe, weights=weights, histograms=False, alpha=False, table=False, rbf=False, simple_grid=True)

	if args.catalogue:
		calibrate(y1v2)

def mkdirs():
	os.system("mkdir -p %s/release/inputs"%config["output_dir"] )
	os.system("mkdir -p %s/release/outputs"%config["output_dir"] )
	os.system("mkdir -p %s/release/rbf"%config["output_dir"] )
	os.system("mkdir -p %s/release/grid"%config["output_dir"] )
	os.system("mkdir -p %s/release/polynomial"%config["output_dir"] )
	os.system("mkdir -p %s/nbc_data"%config["output_dir"] )


def calibrate(data, method="rbf"):

	rbf = (method.lower()=="rbf")

	# Set up the wrappers
	nbc = cal.nbc()
	nbc_disc = cal.nbc()
	nbc_disc.res = data.res
	nbc_bulge = cal.nbc()
	
	# Fit or interpolate
	if method is "polynomial":
		nbc_disc.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc.fits"%config["output_dir"])
		nbc_disc.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc.fits"%config["output_dir"])
		nbc_bulge.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge.fits"%config["output_dir"])
		nbc_bulge.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge.fits"%config["output_dir"])
	elif method is "rbf":
		nbc_disc.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc.fits"%config["output_dir"])
		nbc_bulge.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge.fits"%config["output_dir"])
	elif method is "grid":
		nbc_disc.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-halfcat-disc.fits"%config["output_dir"])[1].read()
		nbc_bulge.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-halfcat-bulge.fits"%config["output_dir"])[1].read()

	# Apply to get a bias correction for each galaxy in the data
	nbc_disc.apply(split_half=0, scheme=method)
	nbc.res = nbc_disc.res
	del(nbc_disc.res)
	print "Done bulge calibration"

	nbc_bulge.res = data.res
	nbc_bulge.apply(split_half=0, scheme=method)
	nbc.res[nbc.res["is_bulge"].astype(bool)] = nbc_bulge.res[nbc.res["is_bulge"].astype(bool)]
	del(nbc_bulge.res)
	print "Done disc calibration"


	if not os.path.exists(os.path.dirname(config["output_dir"])):
		os.system("mkdir -p %s"%os.path.dirname(config["output_dir"]))

	print "Saving calibrated catalogue to %s"%config["output_catalogue"]
	nbc.export(filename="%s/%s"%(config["output_dir"],config["output_catalogue"]))


def diagnostics(y1v2, hoopoe, histograms=True, weights=None, alpha=True, table=True, vssnr=True, vsredshift=True, rbf=True, simple_grid=False):

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
		pl.histograms(["snr", "e", "size", "rgpp", "flux", "psfe", "psf_size"], hoopoe.res, data2=y1v2.res, outdir="%s/release/outputs"%config["output_dir"], weights=weights)

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
	nbc = cal.nbc()
	nbc_disc = cal.nbc()
	nbc_bulge = cal.nbc()
	wt1, wt2 = nbc.get_split_data(hoopoe, weights=weights)

	nbc_disc.load_from_cat(nbc, name="all")
	nbc_bulge.load_from_cat(nbc, name="all")
	

	#### Process disc then bulge runs

	# Fit the binned galaxies. Save one set of calibration data for disc objects, one for bulges
	if table:
		#import pdb ; pdb.set_trace()
		nbc_disc.compute(split_half=1, fit="disc", weights=wt1, rbins=9, sbins=10, rlim=(1.13,3.0), slim=(12,200), table_name="%s/nbc_data/bias_table_hoopoe-v1-halfcat-disc.fits"%config["output_dir"])
		nbc_bulge.compute(split_half=1, fit="bulge", weights=wt1, rbins=9, sbins=10, rlim=(1.13,3.0), slim=(12,200), table_name="%s/nbc_data/bias_table_hoopoe-v1-halfcat-bulge.fits"%config["output_dir"])

		nbc_disc.compute(split_half=0, fit="disc", weights=weights, rbins=9, sbins=10, rlim=(1.13,3.0), slim=(12,200), table_name="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc.fits"%config["output_dir"])
		nbc_bulge.compute(split_half=0, fit="bulge", weights=weights, rbins=9, sbins=10, rlim=(1.13,3.0), slim=(12,200), table_name="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge.fits"%config["output_dir"])

	# Now plot out the points and the resulting smooth fit
	if vssnr:
		if rbf:
			nbc_disc.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc.fits"%config["output_dir"])
			nbc_bulge.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge.fits"%config["output_dir"])
		else:
			nbc_disc.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc.fits"%config["output_dir"])
			nbc_disc.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc.fits"%config["output_dir"])
			nbc_bulge.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge.fits"%config["output_dir"])
			nbc_bulge.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge.fits"%config["output_dir"])

		nbc_disc.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc.fits"%config["output_dir"], bias_name="m", do_half=1, output="%s/m-vs-snr-disc-v1-1.png"%(config["output_dir"]+"/release/"+sub_dir), use_rbf=rbf)
		nbc_disc.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc.fits"%config["output_dir"], bias_name="a", do_half=1, output="%s/alpha-vs-snr-disc-v1-1.png"%(config["output_dir"]+"/release/"+ sub_dir), use_rbf=rbf)
		plt.close()
		nbc_disc.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc.fits"%config["output_dir"], bias_name="m", do_half=2, output="%s/m-vs-snr-disc-v1-2.png"%(config["output_dir"]+"/release/"+ sub_dir), use_rbf=rbf)
		nbc_disc.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-disc.fits"%config["output_dir"], bias_name="a", do_half=2, output="%s/alpha-vs-snr-disc-v1-2.png"%(config["output_dir"]+"/release/"+ sub_dir), use_rbf=rbf)
		plt.close()
		nbc_bulge.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge.fits"%config["output_dir"], bias_name="m", do_half=1, output="%s/m-vs-snr-bulge-v1-2.png"%(config["output_dir"]+"/release/"+ sub_dir), use_rbf=rbf)
		nbc_bulge.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge.fits"%config["output_dir"], bias_name="a", do_half=1, output="%s/alpha-vs-snr-bulge-v1-1.png"%(config["output_dir"]+"/release/"+ sub_dir), use_rbf=rbf)
		plt.close()
		nbc_bulge.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge.fits"%config["output_dir"], bias_name="m", do_half=2, output="%s/m-vs-snr-bulge-v1-1.png"%(config["output_dir"]+"/release/"+ sub_dir), use_rbf=rbf)
		nbc_bulge.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-fullcat-bulge.fits"%config["output_dir"], bias_name="a", do_half=2, output="%s/alpha-vs-snr-bulge-v1-2.png"%(config["output_dir"]+"/release/"+ sub_dir), use_rbf=rbf)
		plt.close()

	if rbf:
		nbc_disc.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-disc.fits"%config["output_dir"])
		nbc_bulge.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-bulge.fits"%config["output_dir"])
	elif simple_grid:
		nbc_disc.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-halfcat-disc.fits"%config["output_dir"])[1].read()
		nbc_bulge.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-halfcat-bulge.fits"%config["output_dir"])[1].read()
	else:
		nbc_disc.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-disc.fits"%config["output_dir"])
		nbc_disc.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-disc.fits"%config["output_dir"])
		nbc_bulge.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-bulge.fits"%config["output_dir"])
		nbc_bulge.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-bulge.fits"%config["output_dir"])

	# Apply to the other half
	if rbf:
		scheme="rbf"
	elif simple_grid:
		scheme="grid"
	else:
		scheme="polynomial"

	nbc_disc.apply(split_half=2, scheme=scheme)
	nbc_bulge.apply(split_half=2, scheme=scheme)

	nbc.get_combined_calibration(nbc_disc,nbc_bulge, split_half=2)

	# Finally save some diagnostic plots in tomographic bins
	if vsredshift: 
		plt.close()
		nbc.redshift_diagnostic(bias="m", label="Uncalibrated", ls="none", nbins=3, fmt=["o","D"], colour="steelblue", weights=wt2)
		nbc.redshift_diagnostic(bias="m", label="Calibrated", ls="none", nbins=3, fmt=["^",">"], apply_calibration=True, colour="purple", weights=wt2)
		plt.ylabel("Multiplicative Bias $m$")
		plt.legend(loc="center right")
		plt.savefig("%s/release/%s/m-vs-redshift-diagnostic-v1.png"%(config["output_dir"],sub_dir))
		plt.close()

		nbc.redshift_diagnostic(bias="alpha", label="Uncalibrated",ls="none", nbins=3, fmt=["o","D"], colour="steelblue", weights=wt2)
		nbc.redshift_diagnostic(bias="alpha", label="Calibrated", ls="none",fmt=["^",">"], nbins=3, apply_calibration=True, colour="purple", weights=wt2)
		plt.ylabel(r"PSF Leakage $\alpha$")
		plt.legend(loc="upper left")
		plt.savefig("%s/release/%s/alpha-vs-redshift-diagnostic-v1.png"%(config["output_dir"],sub_dir))
		plt.close()


if __name__ == "__main__":
	parser = argparse.ArgumentParser(add_help=False)
	parser.add_argument('--config',"-c", type=str, action='store')
	parser.add_argument('--catalogue', action='store_true')
	parser.add_argument('--calculate', action='store_true')

	args = parser.parse_args()

	config = yaml.load(open(args.config))

	mkdirs()
	main(args)

