import numpy as np
import tools.shapes as s
import glob, os, argparse
import fitsio as fi
from tools.im3shape.cosmos import corruption as cpt
import tools.n_of_z as nz

dt = [ ('ra_as', '>f8'), ('dec_as', '>f8'), ('e1', '>f8'), ('e2', '>f8'), ('radius', '>f8'), ('radius_ratio', '>f8'), ('bulge_a', '>f8'), ('disc_a', '>f8'), ('bulge_index', '>f8'), ('disc_index', '>f8'), ('coadd_objects_id', '>i8'), ('time', '>f8'), ('bulge_flux', '>f8'), ('disc_flux', '>f8'), ('flux_ratio', '>f8'), ('snr', '>f8'), ('likelihood', '>f8'), ('n_exposure', '>i8'), ('nparam_varied', '>i8'), ('stamp_size', '>i8'), ('mean_rgpp_rp', '>f8'), ('mean_unmasked_flux_frac', '>f8'), ('mean_hsm_psf_e1_sky', '>f8'), ('mean_hsm_psf_e2_sky', '>f8'), ('mean_hsm_psf_sigma', '>f8'), ('mean_mask_fraction', '>f8'), ('round_snr', '>f8'), ('bands', 'S40'), ('dec', '>f8'), ('tilename', 'S12'), ('ra', '>f8'), ('chi2_pixel', '>f8'), ('modest', '>f8'), ('mag_auto_r', '>f4'), ('flags_r', '>i4'), ('gold_mask', '?'), ('bad_mask', '>f8'), ('error_flag', '>i8'), ('info_flag', '>i8'), ('is_bulge', '>i8'), ('des_bin', '>i8') ]

dt = np.dtype(dt)

# Load the catalogues
def load(filename, header):
	cat = fi.FITS(filename)
	if header!="":
		suffix = "_"+header
	else:
		suffix=""

	hoopoe=s.shapecat()
	hoopoe.res=cat["i3s%s"%suffix].read()
	hoopoe.truth=cat["truth%s"%suffix].read()

	ngal = hoopoe.res.size
	print "Found %3.4fM galaxies."%(ngal/1e6)

	return ngal, hoopoe

def redshift_bins(hoopoe):
    # Assign the COSMOS galaxies to redshift bins
    pz=nz.nofz("/home/samuroff/source.nz")
    pz.generate_interpolators()
    bin_allocation = pz.assign_galaxies_to_bins(hoopoe.truth["cosmos_photoz"])
    return bin_allocation

# Costruct and export the mask array
# Identify the objects with bad COSMOS profiles
def get_mask(filename, hoopoe, ngal):
	print "Constructing object mask"
	gd=cpt.whistleblower()
	select_good = np.invert(np.in1d(hoopoe.truth["cosmos_ident"], gd.blacklisted_ids))

	cut = (hoopoe.res["snr"]>12) & (hoopoe.res["snr"]<200) & (hoopoe.res["mean_rgpp_rp"]>1.13) & (hoopoe.res["mean_rgpp_rp"]<3.00) & (hoopoe.res["info_flag"]==0)

	mask = cut & select_good
	weights = np.ones(ngal)

	if os.path.exists(filename):
		print "file exists %s"%filename
		return 0, 0

	outfits = fi.FITS(filename, "rw")
	outfits.write(mask.astype(int))
	outfits[-1].write_key("EXTNAME", "sel")
	outfits.write(weights)
	outfits[-1].write_key("EXTNAME", "wts")

	outfits.write(hoopoe.res["coadd_objects_id"])
	outfits[-1].write_key("EXTNAME", "coadd_objects_id")

	outfits.close()

	return mask, weights

def catalogue(filename, hoopoe, bin_allocation, ngal):
	out = np.empty(ngal, dtype=dt)
	for col in hoopoe.res.dtype.names:
		if col not in dt.names:
			continue
		else:
			print col
			out[col] = hoopoe.res[col].astype(str(dt[col]))

	out["des_bin"] = bin_allocation

	if os.path.exists(filename):
		print "file exists %s"%filename
		return 0

	outfits = fi.FITS(filename, "rw")
	print "Saving catalogue to %s"%filename

	print "im3shape"
	outfits.write(out)
	outfits[-1].write_key("EXTNAME", "i3s")
	print "truth"
	outfits.write(hoopoe.truth)
	outfits[-1].write_key("EXTNAME", "truth")
	outfits.close()

	print "done"

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--meeting', default="chicago", type=str, action='store')
parser.add_argument('--version', type=int, action='store')
parser.add_argument('--catalogue', default="hoopoe", type=str, action='store')
parser.add_argument('--input', type=str, action='store')
parser.add_argument('--blind', action='store_true')
parser.add_argument('--header', default="", type=str, action='store')

args = parser.parse_args()

if args.blind:
	"BLINDED"
	output = "/share/des/disc8/y1-blinded-combined_cats/"
else:
	"UNBLINDED"
	output = "/share/des/disc8/y1-unblinded-combined_cats/"

os.system("mkdir -p %s"%output)
print output

ngal, hoopoe = load(args.input, args.header)
mask, weight = get_mask("%s/selection_masks/mask-%s_A4_A6-infocuts-%s-i3sv2-catv%.2d-no_Rb_cut.fits"%(output,args.meeting,args.catalogue,args.version), hoopoe, ngal)
bin_allocation = redshift_bins(hoopoe)
catalogue("%s/combined_cats/%s-%s_A4_A6-infocuts-hv2-i3sv2-catv%.2d.fits"%(output,args.meeting,args.catalogue,args.version), hoopoe, bin_allocation, ngal)

print "done all"


