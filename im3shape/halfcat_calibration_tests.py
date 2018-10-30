import numpy as np
import tools.nbc as cal
import tools.shapes as s
import tools.diagnostics as di
import tools.plots as pl
import tools.arrays as arr
import fitsio as fi
import pylab as plt
import os, yaml, argparse, glob, gc

plt.switch_backend("agg")



def do_halfcat_tests(catalogue, weights, params, scheme="grid", split_type="cosmos", sbins=16, rbins=16, smoothing=3):
	# Generate calibrator objects, then split the data
	nbc = cal.nbc()
	nbc_disc = cal.nbc()
	nbc_bulge = cal.nbc()

	wt1, wt2 = nbc.get_split_data(catalogue, weights=weights, method=split_type)
	nbc_disc.load_from_cat(nbc, name="all")
	nbc_bulge.load_from_cat(nbc, name="all")

	gc.collect()
	import pdb ; pdb.set_trace()

	# Do any setting up that needs doing prior to applying the calibration to the catalogue
	nbc_disc, nbc_bulge = setup(nbc_disc, nbc_bulge, split_type, scheme, params["output_dir"], params["sbins"], params["rbins"], smoothing )

	# Apply the calibration from one half of the catalogue to the second part
	nbc_disc.apply(split_half=2, scheme=scheme, names=["m", "a"])
	nbc_bulge.apply(split_half=2, scheme=scheme, names=["m", "a"])
	nbc.combine_bd(nbc_disc,nbc_bulge, split_half=2, names=["m", "c1", "c2"])

	# Finally test for residual bias
	zbins=[ 0.2, 0.43, 0.63, 0.9, 1.3]
	bias0=nbc.redshift_diagnostic(bias="m", label="Uncalibrated", ls="none", nbins=4, fmt=["o","D"], colour="steelblue", weights=wt2, bins=zbins, tophat=False, separate_components=False)
	bias=nbc.redshift_diagnostic(bias="m", label="Calibrated", ls="none", nbins=4, fmt=["^",">"], apply_calibration=True, colour="purple", weights=wt2, bins=zbins, tophat=False, separate_components=False)

	# Save the residual bias numbers to a text file
	np.savetxt("/home/samuroff/shear_pipeline/plot_dump/datavecs/m-vs-z-dvec-%s-halfcat-uncalibrated.txt"%(split_type), np.array(bias0).T, header="z m1 em1 m2 em2 m em")
	np.savetxt("/home/samuroff/shear_pipeline/plot_dump/datavecs/m-vs-z-dvec-%s-halfcat-%s-nbc-sbins%d-rbins%d.txt"%(split_type, scheme, params["sbins"], params["rbins"]), np.array(bias).T, header="z m1 em1 m2 em2 m em")

	print "Done"
	return 0

def setup(nbc_disc, nbc_bulge, split_method, scheme, outdir, sbins, rbins, smoothing):
	"""Read in a grid of bias nodes based on one half of the catalogue.
	   And use it in some way to generate a calibration scheme.
	"""
	if (scheme=="rbf"):
		nbc_disc.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-disc-%dsbins-%drbins.fits"%(outdir, split_method ,sbins, rbins), smoothing=smoothing)
		nbc_bulge.fit_rbf(table="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-bulge-%dsbins-%drbins.fits"%(outdir, split_method ,sbins, rbins), smoothing=smoothing)

	elif (scheme=="grid"):
		nbc_disc.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-disc-%dsbins-%drbins.fits"%(outdir, split_method ,sbins, rbins))[1].read()
		nbc_bulge.bias_grid = fi.FITS("%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-bulge-%dsbins-%drbins.fits"%(outdir, split_method ,sbins, rbins))[1].read()

	elif (scheme=="polynomial"):
		nbc_disc.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-disc-%dsbins-%drbins.fits"%(outdir, split_method ,sbins, rbins))
		nbc_disc.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-disc-%dsbins-%drbins.fits"%(outdir, split_method ,sbins, rbins))
		nbc_bulge.fit("m", table="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-bulge-%dsbins-%drbins.fits"%(outdir, split_method ,sbins, rbins))
		nbc_bulge.fit("a", table="%s/nbc_data/bias_table_hoopoe-v1-%s-halfcat-bulge-%dsbins-%drbins.fits"%(outdir, split_method ,sbins, rbins))

	return nbc_disc, nbc_bulge


