import glob, os, argparse
import fitsio as fi

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--meds', type=str, action='store')
parser.add_argument('--truth', type=str, action='store')
args = parser.parse_args()

truth_dir=args.truth
meds_dir=args.meds
files=glob.glob(meds_dir+"/DES*.fits*")

for i,f in enumerate(files):
	tile=os.path.basename(f)[:12]
	print "%d %s"%(i, tile),
	truth_file=glob.glob("%s/%s*.fz"%(truth_dir,tile))[0]

	tr=fi.FITS(truth_file)["subdetection_objects"].read()
	flux_correction=tr["flux"].sum()

	flux_correction=tr["flux"].sum()/(10000*10000)
	print flux_correction,

	meds=fi.FITS(f, "rw")
	meds["image_cutouts"].write_key("faint_flux_correction", flux_correction)
	meds.close()
	print "done"
