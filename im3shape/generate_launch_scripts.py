import numpy as np
import galsim
import fitsio
import glob, argparse, os
import tools.shapes as s

from string import Template as Tm


def main():
    global args

    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('-f', '--fit', type=str, default="disc", help='fit type (bulge or disc)')
    parser.add_argument('-t', '--template_file', type=str, default="/home/samuroff/shear_pipeline/end-to-end/end-to-end_code/launch/im3shape-template-fornax.sub", help='template for the launch script')
    parser.add_argument('-r', '--repository', type=str, default="/home/samuroff/shear_pipeline", help='code repository')

    args = parser.parse_args()

    os.system("mkdir -p scripts/%s"%args.fit)
    print "Generating launch scripts"

    here = os.getcwd()

    filelist = glob.glob("%s/DES*-meds-*fits.fz"%here)
    print "found %d tiles"%len(filelist)

    for i, f in enumerate(filelist):
    	tile = os.path.basename(f)[:12]
        if "s82" in f:
            continue
        run = f.split("-meds-")[-1].replace(".fz","").replace(".fits","")
        simrun = f.split("-meds-")[0].split("-sim-")[1]
    	print i, tile

        base = open(args.template_file).read()
        base = Tm(base)

        script = base.substitute(TILE=tile, GROUP=run, SIMRUN=simrun, DIRNAME=here, FIT=args.fit, SHEARPIPE=args.repository, PBS_ARRAYID="$PBS_ARRAYID",PBS_JOBID="$PBS_JOBID")
        if ".fz" not in f:
            script = script.replace(".fz","")

        out="%s/scripts/%s"%(here,args.fit)

        outfile = open("%s/%s.sub"%(out,tile), "w")
        outfile.write(script)
        outfile.close()

    print "Done"


main()
