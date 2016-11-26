import tools.nbc as cal
import tools.shapes as s
import tools.diagnostics as di
import tools.plots as pl
import fitsio as fi
import pylab as plt
import os


# Input im3shape results on data
i3s_dir = "/share/des/disc7/samuroff/des/im3shape-v2-infocuts.fits" 
# Where to save the output plots and bias table
plots_dir = "/home/samuroff/shear_pipeline/bias_calibration/v1.0/release/"
table_dir ="/home/samuroff/shear_pipeline/bias_calibration/v1.0/nbc_data/"
# The final calibrated catalogue -  the end result should just be the i3s_dir file saved with some extra columns
output_catalogue = "/home/samuroff/shear_pipeline/bias_calibration/v1.0/y1a1-im3shape-r-1-2-2.fits"

plt.switch_backend("agg")

def main():
	im3shape_columns = ["e1", "e2", "mean_psf_e1_sky", "mean_psf_e2_sky", "snr", "mean_rgpp_rp", "radius", "coadd_objects_id", "mean_flux", "n_exposure", "stamp_size", "info_flag", "is_bulge", "tilename"]
	truth_columns = ['DES_id', 'cosmos_ident', 'cosmos_photoz', 'sextractor_pixel_offset', 'true_g1', 'true_g2', 'intrinsic_e1', 'intrinsic_e2', 'ra', 'dec', 'hlr', 'mag', 'flux']
	# Load the y1 data
	y1v2 = s.shapecat(res=i3s_dir)
	y1v2.load(truth=False, prune=True, cols=[im3shape_columns,truth_columns])
	y1v2.res=y1v2.res[y1v2.res["info_flag"]==0]

	# And the simulation results
	hoopoe = s.shapecat(res="/share/des/disc6/samuroff/y1/hoopoe/y1_collated/results/disc-fits/main",truth="/share/des/disc6/samuroff/y1/hoopoe/y1_collated/truth")
	hoopoe.load(truth=True, cols=[im3shape_columns,truth_columns])
	
	diagnostics(y1v2, hoopoe, rbf=True)
	diagnostics(y1v2, hoopoe, histograms=False, alpha=False, table=False, rbf=False)
	calibrate(y1v2, hoopoe)

def calibrate(data, hoopoe, method="rbf"):

	rbf = (method.lower()=="rbf")

	# Set up the wrappers
	nbc = cal.nbc()
	nbc_disc = cal.nbc()
	nbc_bulge = cal.nbc()
	nbc.get_split_data(data)
	nbc_disc.load_from_cat(nbc, name="res")
	nbc_bulge.load_from_cat(nbc, name="res")

	# Evaluate at fixed points
	nbc_disc.compute(split_half=0, fit="disc", slim=(10,350), table_name="%s/bias_table_hoopoe-v1-fullcat-disc.fits"%table_dir)
	nbc_bulge.compute(split_half=0, fit="bulge", slim=(10,350), table_name="%s/bias_table_hoopoe-v1-fullcat-bulge.fits"%table_dir)

	# Fit or interpolate
	if method is "polynomial":
		nbc_disc.fit("m", table="%s/bias_table_hoopoe-v1-fullcat-disc.fits"%table_dir)
		nbc_disc.fit("a", table="%s/bias_table_hoopoe-v1-fullcat-disc.fits"%table_dir)
		nbc_bulge.fit("m", table="%s/bias_table_hoopoe-v1-fullcat-bulge.fits"%table_dir)
		nbc_bulge.fit("a", table="%s/bias_table_hoopoe-v1-fullcat-bulge.fits"%table_dir)
	elif method is "rbf":
		nbc_disc.fit_rbf(table="%s/bias_table_hoopoe-v1-halfcat-disc.fits"%table_dir)
		nbc_bulge.fit_rbf(table="%s/bias_table_hoopoe-v1-halfcat-bulge.fits"%table_dir)

	# Apply to get a bias correction for eacg galaxy in the data
	nbc_disc.apply(split_half=0, use_rbf=rbf)
	nbc_bulge.apply(split_half=0, use_rbf=rbf)
	nbc.get_combined_calibration(nbc_disc,nbc_bulge, split_half=0)

	if not os.path.exists(os.path.dirname(output_catalogue)):
		os.system("mkdir -p %s"%os.path.dirname(output_catalogue))

	nbc.export(filename=output_catalogue)


def diagnostics(y1v2, hoopoe, histograms=True, alpha=True, table=True, vssnr=True, vsredshift=True, rbf=True):

	if not os.path.exists(plots_dir):
		print "Warning - output plots directory does not exist"

	# First make some diagnostic distributions of the observable parameters relevant for the NBC
	if histograms:
		pl.histograms_vs_input(["e", "size", "flux"], hoopoe.truth, data2=y1v2.res, outdir="%s/inputs"%plots_dir)
		pl.histograms(["snr", "e", "size", "rgpp", "flux", "psfe", "psf_size"], hoopoe.res, data2=y1v2.res, outdir="%s/outputs"%plots_dir)

	if alpha:
		# Global fit of alpha to compare with previous runs
		plt.close()
		b = di.get_alpha(hoopoe.res, hoopoe.res, nbins=20, xlim=(-0.015,0.025), binning="equal_number", use_weights=False, visual=True)
		plt.subplots_adjust(wspace=0, hspace=0, top=0.95, bottom=0.1)
		plt.title("Y1 v1 Sim, PSF v01")

		os.system("mkdir -p %s/alpha"%plots_dir)
		plt.savefig("%s/alpha/alphaplot-e-vs-epsf-linfit.png"%plots_dir)
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
		nbc_disc.compute(split_half=1, fit="disc", rbins=10, sbins=10, rlim=(1.0,2.0), slim=(10,200), table_name="%s/bias_table_hoopoe-v1-halfcat-disc.fits"%table_dir)
		nbc_bulge.compute(split_half=1, fit="bulge", rbins=10, sbins=10, rlim=(1.0,2.0), slim=(10,200), table_name="%s/bias_table_hoopoe-v1-halfcat-bulge.fits"%table_dir)
	
	if rbf:
		nbc_disc.fit_rbf(table="%s/bias_table_hoopoe-v1-halfcat-disc.fits"%table_dir)
		nbc_bulge.fit_rbf(table="%s/bias_table_hoopoe-v1-halfcat-bulge.fits"%table_dir)
	else:
		nbc_disc.fit("m", table="%s/bias_table_hoopoe-v1-halfcat-disc.fits"%table_dir)
		nbc_disc.fit("a", table="%s/bias_table_hoopoe-v1-halfcat-disc.fits"%table_dir)
		nbc_bulge.fit("m", table="%s/bias_table_hoopoe-v1-halfcat-bulge.fits"%table_dir)
		nbc_bulge.fit("a", table="%s/bias_table_hoopoe-v1-halfcat-bulge.fits"%table_dir)
	

	# Now plot out the points and the resulting smooth fit
	if vssnr:
		nbc_disc.bias_fit_vs_pts(table="%s/bias_table_hoopoe-v1-halfcat-disc.fits"%table_dir, bias_name="m", do_half=1, output="%s/m-vs-snr-disc-v1-1.png"%(plots_dir+"/rbf"*rbf), use_rbf=rbf)
		nbc_disc.bias_fit_vs_pts(table="%s/bias_table_hoopoe-v1-halfcat-disc.fits"%table_dir, bias_name="a", do_half=1, output="%s/alpha-vs-snr-disc-v1-1.png"%(plots_dir+"/rbf"*rbf), use_rbf=rbf)
		plt.close()
		nbc_disc.bias_fit_vs_pts(table="%s/bias_table_hoopoe-v1-halfcat-disc.fits"%table_dir, bias_name="m", do_half=2, output="%s/m-vs-snr-disc-v1-2.png"%(plots_dir+"/rbf"*rbf), use_rbf=rbf)
		nbc_disc.bias_fit_vs_pts(table="%s/bias_table_hoopoe-v1-halfcat-disc.fits"%table_dir, bias_name="a", do_half=2, output="%s/alpha-vs-snr-disc-v1-2.png"%(plots_dir+"/rbf"*rbf), use_rbf=rbf)
		plt.close()
		nbc_bulge.bias_fit_vs_pts(table="%s/bias_table_hoopoe-v1-halfcat-bulge.fits"%table_dir, bias_name="m", do_half=1, output="%s/m-vs-snr-bulge-v1-2.png"%(plots_dir+"/rbf"*rbf), use_rbf=rbf)
		nbc_bulge.bias_fit_vs_pts(table="%s/bias_table_hoopoe-v1-halfcat-bulge.fits"%table_dir, bias_name="a", do_half=1, output="%s/alpha-vs-snr-bulge-v1-1.png"%(plots_dir+"/rbf"*rbf), use_rbf=rbf)
		plt.close()
		nbc_bulge.bias_fit_vs_pts(table="%s/bias_table_hoopoe-v1-halfcat-bulge.fits"%table_dir, bias_name="m", do_half=2, output="%s/m-vs-snr-bulge-v1-1.png"%(plots_dir+"/rbf"*rbf), use_rbf=rbf)
		nbc_bulge.bias_fit_vs_pts(table="%s/bias_table_hoopoe-v1-halfcat-bulge.fits"%table_dir, bias_name="a", do_half=2, output="%s/alpha-vs-snr-bulge-v1-2.png"%(plots_dir+"/rbf"*rbf), use_rbf=rbf)
		plt.close()

	# Apply to the other half
	nbc_disc.apply(split_half=2, use_rbf=rbf)
	nbc_bulge.apply(split_half=2, use_rbf=rbf)

	nbc.get_combined_calibration(nbc_disc,nbc_bulge, split_half=2)

	# Finallt save some diagnostic plots in tomographic bins
	if vsredshift: 
		nbc.redshift_diagnostic(bias="m", label="Uncalibrated", ls="none", nbins=3, fmt=["o","D"], colour="steelblue")
		nbc.redshift_diagnostic(bias="m", label="Calibrated", ls="none", nbins=3, fmt=["^",">"], apply_calibration=True, colour="purple")
		plt.ylabel("Multiplicative Bias $m$")
		plt.legend(loc="lower left")
		plt.savefig("%s/m-vs-redshift-diagnostic-v1.png"%(plots_dir+rbf*"/rbf"))
		plt.close()

		nbc.redshift_diagnostic(bias="alpha", label="Uncalibrated",ls="none", nbins=3, fmt=["o","D"], colour="steelblue")
		nbc.redshift_diagnostic(bias="alpha", label="Calibrated", ls="none",fmt=["^",">"], nbins=3, apply_calibration=True, colour="purple")
		plt.ylabel(r"PSF Leakage $\alpha$")
		plt.legend(loc="lower right")
		plt.savefig("%s/alpha-vs-redshift-diagnostic-v1.png"%(plots_dir+rbf*"/rbf"))
		plt.close()










#main()

