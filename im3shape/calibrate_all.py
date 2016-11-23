import tools.nbc as cal
import tools.shapes as s
import tools.diagnostics as di
import tools.plots as pl
import fitsio as fi
import pylab as plt
import os

rbf=False

i3s_dir = "/share/des/disc2/y1/y1a1-im3shape-r-1-1-1.fits"
#"/share/des/disc7/samuroff/des/im3shape-v2-infocuts.fits" 
plots_dir = "/home/samuroff/shear_pipeline/bias_calibration/v2.2/release/"
table_dir ="/home/samuroff/shear_pipeline/bias_calibration/v2.2/nbc_data/"

plt.switch_backend("agg")

def main():
	# Load the y1 data
	y1v2 = s.shapecat(res=i3s_dir)
	y1v2.load(truth=False)
	y1v2.res=y1v2.res[y1v2.res["info_flag"]==0]

	hoopoe = s.shapecat(res="/share/des/disc6/samuroff/y1/hoopoe/y1_collated/results/disc-fits/main",truth="/share/des/disc6/samuroff/y1/hoopoe/y1_collated/truth")
	#hoopoe = s.shapecat(res="/share/des/disc6/samuroff/y1/hoopoe/y1a1-v2.2_10/results/full/results/blinded/bord/main/",truth="/share/des/disc3/samuroff/y1/sims/v2.2/y1a1_16tiles/truth/")
	hoopoe.load(truth=True)
	
	diagnostics(y1v2,hoopoe)
	#calibration()

def diagnostics(y1v2,hoopoe):

	# First make some diagnostic distributions of the observable parameters relevant for the NBC
	pl.histograms_vs_input(["redshift", "e", "size", "flux"], hoopoe.truth, data2=y1v2.res, outdir="%s/inputs"%plots_dir)
	pl.histograms(["snr", "e", "size", "rgpp", "flux", "psfe", "psf_size"], hoopoe.res, data2=y1v2.res, outdir="%s/outputs"%plots_dir)

	# Global fit of alpha to compare with previous runs
	plt.close()
	b = di.get_alpha(hoopoe.res, nbins=20, xlim=(-0.025,0.025), binning="equal_number", use_weights=False, visual=True)
	os.system("mkdir -p %s/alpha"%plots_dir)
	plt.savefig("%s/alpha/alphaplot-e-vs-epsf-linfit.png"%plots_dir)
	plt.close()

	# Now compute the calibration using a random half of the simulation
	nbc = cal.nbc()
	nbc.get_split_data(hoopoe)
	nbc_disc = cal.nbc()
	nbc_disc.get_split_data(hoopoe)
	nbc_bulge = cal.nbc()
	nbc_bulge.get_split_data(hoopoe)

	# Fit the disc galaxies and save one set of calibration data
	#nbc_disc.compute(split_half=1, fit="bord", slim=(10,350), table_name="%s/bias_table_hoopoe-v2-halfcat-disc.fits"%table_dir)
	
	if rbf:
		nbc_disc.fit_rbf(table="%s/bias_table_hoopoe-v2-halfcat-disc.fits"%table_dir)
	else:
		nbc_disc.fit("m", table="%s/bias_table_hoopoe-v2-halfcat-disc.fits"%table_dir)
		nbc_disc.fit("a", table="%s/bias_table_hoopoe-v2-halfcat-disc.fits"%table_dir)

	nbc_disc.apply(split_half=2, use_rbf=rbf)

	# Fit the bulge galaxies and save another set of calibration data
	#nbc_bulge.compute(split_half=1, fit="bulge", table_name="%s/bias_table_hoopoe-v2-halfcat-bulge.fits"%table_dir)
	#nbc_bulge.fit("m", table_name="%s/bias_table_hoopoe-v2-halfcat-bulge.fits"%table_dir)
	#nbc_bulge.fit("a", table_name="%s/bias_table_hoopoe-v2-halfcat-bulge.fits"%table_dir)
	#nbc_bulge.apply(split_half=2)
#
	#nbc.get_combined_calibration(nbc_disc,nbc_bulge)

	# Now plot out the fit points and the resulting polynomial
	nbc_disc.bias_fit_vs_pts(table="%s/bias_table_hoopoe-v2-halfcat-disc.fits"%table_dir, bias_name="m", output="%s/m-vs-snr-disc-v2.png"%(plots_dir+"/rbf"*rbf), use_rbf=rbf)
	nbc_disc.bias_fit_vs_pts(table="%s/bias_table_hoopoe-v2-halfcat-disc.fits"%table_dir, bias_name="a", output="%s/alpha-vs-snr-disc-v2.png"%(plots_dir+"/rbf"*rbf), use_rbf=rbf)
	plt.close()
	nbc=nbc_disc

	# Save some diagnostic plots in tomographic bins
	nbc.redshift_diagnostic(bias="m", label="Uncalibrated", nbins=3, fmt="D", colour="steelblue")
	nbc.redshift_diagnostic(bias="m", label="Calibrated", nbins=3, apply_calibration=True, colour="purple")
	plt.ylabel("Multiplicative Bias $m$")
	plt.legend(loc="lower left")
	plt.savefig("%s/m-vs-redshift-diagnostic-v2.png"%(plots_dir+rbf*"/rbf"))
	plt.close()

	nbc.redshift_diagnostic(bias="alpha", label="Uncalibrated", nbins=3, fmt="D", colour="steelblue")
	nbc.redshift_diagnostic(bias="alpha", label="Calibrated", nbins=3, apply_calibration=True, colour="purple")
	plt.ylabel(r"PSF Leakage $\alpha$")
	plt.legend(loc="lower right")
	plt.savefig("%s/alpha-vs-redshift-diagnostic-v2.png"%(plots_dir+rbf*"/rbf"))
	plt.close()










main()

