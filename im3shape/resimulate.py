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
    parser.add_argument('-f', '--field', type=str, default="*", help='Simulation run to resimulate.')
    parser.add_argument('--mpi',  action='store_true', help='Split the tiles by MPI rank' )
    args = parser.parse_args()

    if args.mpi:
        print "Setting up MPI."
        import mpi4py.MPI
        rank = mpi4py.MPI.COMM_WORLD.Get_rank()
        size = mpi4py.MPI.COMM_WORLD.Get_size()
    else:
        print "Not using MPI (set the --mpi flag if you do want to parallelise the calculation)"
        rank = 0
        size = 1

    os.system("mkdir -p %s"%args.output)

    print "Will remove %s bias"%args.mode
    print "Resimulated tiles will be written to %s"%args.output

    filelist = glob.glob("%s/DES*%s*fits*.fz"%(args.input,args.field))
    print "found %d tiles"%len(filelist)

    for i, f in enumerate(filelist):
    	if os.path.exists("%s/%s"%(args.output,os.path.basename(f))):
            print "file exists."
            continue
        if i%size!=rank:
            continue
    	print i, f

    	new_meds = f
        #"%s/%s"%(args.output,os.path.basename(f))

    	truth_file = os.path.basename(f).replace("-meds-","-truth-")

    	cat = s.shapecat(res=None, truth="%s/%s"%(args.truth, truth_file))
    	if args.mode=="model":
            load_truth = True
            cat.load(res=False, truth=True)

    	meds = s.meds_wrapper(new_meds, update=False)

    	if (args.mode=="model"):
            try:
                meds.remove_model_bias(cat, silent=True, outdir=args.output, noise=False, neighbours=False)
            except:
                continue
    	elif (args.mode=="neighbour"):
    		meds.remove_neighbours(silent=False, outdir=args.output, noise=True)
    	elif (args.mode=="noise"):
            meds.remove_noise(silent=True, outdir=args.output)
        elif (args.mode=="neighbour_noise"):
            meds.remove_neighbours(silent=False, outdir=args.output, noise=False)


print "Done"


main()


def reimpose_neighbour_masks(source_galaxies, source_masks, target, compressed=True):
    """Direct this function towards a set of neighbour free MEDS files"""
    if compressed:
        suffix = ".fz"
    else:
        suffix = ""

    files = glob.glob("%s/DES*.fits%s"%(source, suffix))

    print "Will remove neighbour masking in %d MEDS files."

    for i, f in files:
        tile = os.basename(f)[:12]
        print i, tile
        m = s.meds_wrapper(f)



