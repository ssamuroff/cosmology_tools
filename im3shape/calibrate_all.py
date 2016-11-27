import numpy as np
import tools.nbc as cal
import tools.shapes as s
import tools.diagnostics as di
import tools.plots as pl
import fitsio as fi
import pylab as plt
import os, yaml, argparse

plt.switch_backend("agg")

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--config',"-c", type=str, action='store')
parser.add_argument('--catalogue', action='store_true')
parser.add_argument('--calculate', action='store_true')

args = parser.parse_args()

config = yaml.load(open(args.config))

def main(args):
	im3shape_columns = ["e1", "e2", "mean_psf_e1_sky", "mean_psf_e2_sky", "mean_psf_fwhm", "snr", "mean_rgpp_rp", "radius", "coadd_objects_id", "mean_flux", "n_exposure", "stamp_size", "info_flag", "is_bulge", "tilename"]
	truth_columns = ['DES_id', 'cosmos_ident', 'cosmos_photoz', 'sextractor_pixel_offset', 'true_g1', 'true_g2', 'intrinsic_e1', 'intrinsic_e2', 'ra', 'dec', 'hlr', 'mag', 'flux']

	# Load the y1 data
	if args.calculate or args.catalogue:
		y1v2 = s.shapecat(res=config["i3s_dir"])
		y1v2.load(truth=False, prune=True, cols=[im3shape_columns,truth_columns])
		y1v2.res=y1v2.res[y1v2.res["info_flag"]==0]
	

	# And the simulation results
	if args.calculate:
		hoopoe = s.shapecat(res="%s/bord-fits/main"%config["hoopoe_dir"] ,truth="%s/truth"%config["hoopoe_dir"])
		hoopoe.load(truth=True, cols=[im3shape_columns,truth_columns])
		#hoopoe.res=fi.FITS("/share/des/disc7/samuroff/des/hoopoe-y1a1-v02-bord.fits")["i3s"].read()
		#hoopoe.truth=fi.FITS("/share/des/disc7/samuroff/des/hoopoe-y1a1-v02-bord.fits")["truth"].read()
		sel = np.isfinite(hoopoe.res["mean_psf_e1_sky"]) & np.isfinite(hoopoe.res["mean_psf_e2_sky"])
		hoopoe.res = hoopoe.res[sel]
		hoopoe.truth = hoopoe.truth[sel]

		diagnostics(y1v2, hoopoe, vssnr=config["output"]["snr"], vsredshift=config["output"]["redshift"], table=config["output"]["tables"], alpha=config["output"]["alpha"], histograms=config["output"]["histograms"], rbf=True)
		diagnostics(y1v2, hoopoe, histograms=False, alpha=False, table=False, rbf=False)

	if args.catalogue:
		calibrate(y1v2)

def mkdirs():
	os.system("mkdir -p %s/release/inputs"%config["output_dir"] )
	os.system("mkdir -p %s/release/outputs"%config["output_dir"] )
	os.system("mkdir -p %s/release/rbf"%config["output_dir"] )
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
		nbc_disc.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-disc.fits"%config["output_dir"])
		nbc_bulge.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-bulge.fits"%config["output_dir"])

	# Apply to get a bias correction for each galaxy in the data
	nbc_disc.apply(split_half=0, use_rbf=rbf)
	nbc.res = nbc_disc.res
	del(nbc_disc.res)
	print "Done bulge calibration"

	nbc_bulge.res = data.res
	nbc_bulge.apply(split_half=0, use_rbf=rbf)
	nbc.res[nbc.res["is_bulge"].astype(bool)] = nbc_bulge.res[nbc.res["is_bulge"].astype(bool)]
	del(nbc_bulge.res)
	print "Done disc calibration"


	if not os.path.exists(os.path.dirname(config["output_dir"])):
		os.system("mkdir -p %s"%os.path.dirname(config["output_dir"]))

	print "Saving calibrated catalogue to %s"%config["output_catalogue"]
	nbc.export(filename="%s/%s"%(config["output_dir"],config["output_catalogue"]))


def diagnostics(y1v2, hoopoe, histograms=True, alpha=True, table=True, vssnr=True, vsredshift=True, rbf=True):

	if not os.path.exists(config["output_dir"]):
		print "Warning - output plots directory does not exist"

	# First make some diagnostic distributions of the observable parameters relevant for the NBC
	if histograms:
		pl.histograms_vs_input(["e", "size", "flux"], hoopoe.truth, data2=y1v2.res, outdir="%s/release/inputs"%config["output_dir"])
		pl.histograms(["snr", "e", "size", "rgpp", "flux", "psfe", "psf_size"], hoopoe.res, data2=y1v2.res, outdir="%s/release/outputs"%config["output_dir"])

	if alpha:
		# Global fit of alpha to compare with previous runs
		plt.close()
		b = di.get_alpha(hoopoe.res, hoopoe.res, nbins=20, xlim=(-0.015,0.025), binning="equal_number", use_weights=False, visual=True)
		plt.subplots_adjust(wspace=0, hspace=0, top=0.95, bottom=0.1)
		plt.title("Y1 v2 Sim, PSF v02")

		os.system("mkdir -p %s/release/alpha"%config["output_dir"])
		plt.savefig("%s/release/alpha/alphaplot-e-vs-epsf-linfit.png"%config["output_dir"])
		plt.close()


	#Exit here if no calibration diagnostics are needed
	if (not table) and (not vssnr) and (not vsredshift):
		return 0

	# Now compute the calibration using a random half of the simulation
	nbc = cal.nbc()
	nbc_disc = cal.nbc()
	nbc_bulge = cal.nbc()
	nbc.get_split_data(hoopoe)

	nbc_disc.load_from_cat(nbc, name="all")
	nbc_bulge.load_from_cat(nbc, name="all")
	

	#### Process disc then bulge runs

	# Fit the binned galaxies. Save one set of calibration data for disc objects, one for bulges
	if table:
		nbc_disc.compute(split_half=1, fit="disc", rbins=10, sbins=10, rlim=(0.9,4.5), slim=(10,200), table_name="%s/nbc_data/bias_table_hoopoe-v1-halfcat-disc.fits"%config["output_dir"])
		nbc_bulge.compute(split_half=1, fit="bulge", rbins=10, sbins=10, rlim=(0.9,4.5), slim=(10,200), table_name="%s/nbc_data/bias_table_hoopoe-v1-halfcat-bulge.fits"%config["output_dir"])
	
	if rbf:
		nbc_disc.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-disc.fits"%config["output_dir"])
		nbc_bulge.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-bulge.fits"%config["output_dir"])
	else:
		nbc_disc.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-disc.fits"%config["output_dir"])
		nbc_disc.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-disc.fits"%config["output_dir"])
		nbc_bulge.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-bulge.fits"%config["output_dir"])
		nbc_bulge.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-bulge.fits"%config["output_dir"])
	

	# Now plot out the points and the resulting smooth fit
	if vssnr:
		nbc_disc.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-disc.fits"%config["output_dir"], bias_name="m", do_half=1, output="%s/m-vs-snr-disc-v1-1.png"%(config["output_dir"]+"/release"+"/rbf"*rbf), use_rbf=rbf)
		nbc_disc.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-disc.fits"%config["output_dir"], bias_name="a", do_half=1, output="%s/alpha-vs-snr-disc-v1-1.png"%(config["output_dir"]+"/release"+"/rbf"*rbf), use_rbf=rbf)
		plt.close()
		nbc_disc.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-disc.fits"%config["output_dir"], bias_name="m", do_half=2, output="%s/m-vs-snr-disc-v1-2.png"%(config["output_dir"]+"/release"+"/rbf"*rbf), use_rbf=rbf)
		nbc_disc.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-disc.fits"%config["output_dir"], bias_name="a", do_half=2, output="%s/alpha-vs-snr-disc-v1-2.png"%(config["output_dir"]+"release"+"/rbf"*rbf), use_rbf=rbf)
		plt.close()
		nbc_bulge.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-bulge.fits"%config["output_dir"], bias_name="m", do_half=1, output="%s/m-vs-snr-bulge-v1-2.png"%(config["output_dir"]+"/release"+"/rbf"*rbf), use_rbf=rbf)
		nbc_bulge.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-bulge.fits"%config["output_dir"], bias_name="a", do_half=1, output="%s/alpha-vs-snr-bulge-v1-1.png"%(config["output_dir"]+"/release"+"/rbf"*rbf), use_rbf=rbf)
		plt.close()
		nbc_bulge.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-bulge.fits"%config["output_dir"], bias_name="m", do_half=2, output="%s/m-vs-snr-bulge-v1-1.png"%(config["output_dir"]+"/release"+"/rbf"*rbf), use_rbf=rbf)
		nbc_bulge.bias_fit_vs_pts(table="%s/nbc_data/bias_table_hoopoe-v1-halfcat-bulge.fits"%config["output_dir"], bias_name="a", do_half=2, output="%s/alpha-vs-snr-bulge-v1-2.png"%(config["output_dir"]+"/release"+"/rbf"*rbf), use_rbf=rbf)
		plt.close()

	# Apply to the other half
	nbc_disc.apply(split_half=2, use_rbf=rbf)
	nbc_bulge.apply(split_half=2, use_rbf=rbf)

	nbc.get_combined_calibration(nbc_disc,nbc_bulge, split_half=2)

	# Finally save some diagnostic plots in tomographic bins
	if vsredshift: 
		nbc.redshift_diagnostic(bias="m", label="Uncalibrated", ls="none", nbins=3, fmt=["o","D"], colour="steelblue")
		nbc.redshift_diagnostic(bias="m", label="Calibrated", ls="none", nbins=3, fmt=["^",">"], apply_calibration=True, colour="purple")
		plt.ylabel("Multiplicative Bias $m$")
		plt.legend(loc="lower left")
		plt.savefig("%s/release/%s/m-vs-redshift-diagnostic-v1.png"%(config["output_dir"],rbf*"/rbf"))
		plt.close()

		nbc.redshift_diagnostic(bias="alpha", label="Uncalibrated",ls="none", nbins=3, fmt=["o","D"], colour="steelblue")
		nbc.redshift_diagnostic(bias="alpha", label="Calibrated", ls="none",fmt=["^",">"], nbins=3, apply_calibration=True, colour="purple")
		plt.ylabel(r"PSF Leakage $\alpha$")
		plt.legend(loc="lower right")
		plt.savefig("%s/release/%s/alpha-vs-redshift-diagnostic-v1.png"%(config["output_dir"],rbf*"/rbf"))
		plt.close()









mkdirs()
main(args)

