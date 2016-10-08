import numpy as np
import galsim
import fitsio
import glob, argparse, os
import tools.shapes as s


def main():
    global args

    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('-v', '--verbosity', type=int, action='store', default=2, choices=(0, 1, 2, 3), help='integer verbosity level: min=0, max=3 [default=2]')
    parser.add_argument('-i', '--input', type=str, default=".", help='Directory to source MEDS files')

    args = parser.parse_args()

    filelist = glob.glob("%s/DES*-meds-*fits*.fz"%args.input)
    print "found %d tiles"%len(filelist)

    for i, f in enumerate(filelist):
    	#os.system("cp -rf %s %s"%(f, args.output))
    	if args.verbosity>0:
            print i, f

    	meds = s.meds_wrapper(f, update=False)
        meds.download_source_images()


print "Done"


main()
