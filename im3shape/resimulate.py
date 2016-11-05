import numpy as np
import galsim
import fitsio
import glob, argparse, os
import tools.shapes as s


def main():
    global args

    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-o', '--output', type=str, default="resimulated", help='Directory to write to.')
    parser.add_argument('-i', '--input', type=str, default=".", help='Directory to source MEDS files')
    parser.add_argument('-t', '--truth', type=str, default="truth", help='Directory to source truth tables')
    parser.add_argument('-m', '--mode', type=str, default="model", help='type of bias to remove.')
    args = parser.parse_args()

    os.system("mkdir -p %s"%args.output)
    print "Resimulated tiles will be written to %s"%args.output

    filelist = glob.glob("%s/DES*-meds-*fits*.fz"%args.input)
    print "found %d tiles"%len(filelist)

    for i, f in enumerate(filelist):
    	if os.path.exists("%s/%s"%(args.output,f)):
            print "file exists."
            continue
    	print i, f

    	new_meds = f
        #"%s/%s"%(args.output,os.path.basename(f))

    	truth_file = os.path.basename(f).replace("-meds-","-truth-")

    	cat = s.shapecat(res=None, truth="%s/%s"%(args.truth, truth_file))
    	cat.load(res=False, truth=True)

    	meds = s.meds_wrapper(new_meds, update=False)

    	if (args.mode=="model"):
            meds.remove_model_bias(cat, silent=True, outdir=args.output)
    	elif (args.mode=="neighbour"):
    		meds.remove_neighbours(silent=False, outdir=args.output, noise=True)
    	elif (args.mode=="noise"):
            meds.remove_noise(silent=True, outdir=args.output)


print "Done"


main()


def remove_neighbour_masks(source, target, compressed=True):
    if compressed:
        suffix = ".fz"
    else:
        suffix = ""

    files = glob.glob("%s/DES*.fits%s"%(source, suffix))

    for i, f in files:
        tile = os.basename(f)[:12]
        print i, tile
        m = s.meds_wrapper(f)



