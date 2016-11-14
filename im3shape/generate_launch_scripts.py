import numpy as np
import galsim
import astropy.io.fits as fitsio
import glob, argparse, os
#import tools.shapes as s

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


#main()

class meds_list:
    def __init__(self, results="r/disc/", meds="."):
        if meds==".":
            meds = os.getcwd() 
        print "Initialised"
        print "existing results: %s"%results
        print "MEDS files: %s"%meds
        self.meds_dir = meds
        self.results_dir = results

    def get_tiles_to_do(self):
        
        meds_tiles = glob.glob("%s/DES*.fz"%self.meds_dir)
        results = glob.glob("%s/DES*/1.main.txt"%self.results_dir)
        results_tiles = [os.path.dirname(r)[-12:] for r in results if (os.path.getsize(r)>0)]
        print "%d MEDS files"%len(meds_tiles)
        print "%d tiles with im3shape results"%len(results)

        self.todo = []

        for i, f in enumerate(meds_tiles):
            tile = os.path.basename(f)[:12]
            print i, tile,

            if (tile in results_tiles):
                print "done"
                continue
            else:
                print "to do"
                self.todo.append(f)

        print "Found %d tiles to process"%len(self.todo)

    def write(self, filename="meds_list.ini"):
        print "writing to list %s"%filename
        out = open(filename, "wa")
        for file in self.todo:
            out.write("%s \n"%file)
        out.close()
