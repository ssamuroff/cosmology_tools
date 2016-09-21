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
    args = parser.parse_args()

    os.system("mkdir -p %s"%args.output)
    print "Resimulated tiles will be written to %s"%args.output

    filelist = glob.glob("%s/DES*-meds-*fits*"%args.input)
    print "found %d tiles"%len(filelist)

    for i, f in enumerate(filelist):
    	os.system("cp %s %s"%(f, args.output))
    	print i, f

    	new_meds = "%s/%s"%(args.output,os.path.basename(f))

    	truth_file = os.path.basename(f).replace("-meds-","-truth-")

    	cat = s.shapecat(res=None, truth="%s/%s"%(args.truth, truth_file))
    	cat.load(res=False, truth=True)

    	meds = s.meds_wrapper(new_meds, update=True)
    	meds.remove_model_bias(cat)

print "Done"


main()
