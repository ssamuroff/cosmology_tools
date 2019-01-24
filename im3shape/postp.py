import glob
import argparse
import re
import errno
import sys

class FakeMPI(object):
    def Get_rank(self):
        return 0
    def Get_size(self):
        return 1
#world = FakeMPI()
from mpi4py.MPI import COMM_WORLD as world
report=False

from des_post.postprocess import process_text

tile_pattern = re.compile(r'DES[0-9][0-9][0-9][0-9][+-][0-9][0-9][0-9][0-9]')
tile_band_pattern = re.compile(r'DES[0-9][0-9][0-9][0-9][+-][0-9][0-9][0-9][0-9][_-][ugrizy]')
run_pattern = re.compile(r'[0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]_DES[0-9][0-9][0-9][0-9][+-][0-9][0-9][0-9][0-9]')
great_des_pattern=re.compile(r"nbc(.)+\.meds\.([0-9][0-9][0-9])\.g([0-9][0-9])\.fits")


def find_tilename(name):
    m = tile_pattern.search(name)
    if m is None:
        raise ValueError("Unknown tile {}".format(name))
    return m.group()


import os

def process_dir(indir, outdir):
    main_text_files = glob.glob("{0}/main/*.txt".format(indir))
    rank = world.Get_rank()
    size = world.Get_size()
    main_text_files_2 = []
    for m in main_text_files:
        tilename = find_tilename(m)
        out_main = "{0}/main/{1}.fits".format(outdir, tilename)
        out_epoch = "{0}/epoch/{1}.fits".format(outdir, tilename)
        if not (os.path.exists(out_main) and os.path.exists(out_epoch)):
            main_text_files_2.append(m)
    main_text_files=main_text_files_2
    print "{0} files left to do".format(len(main_text_files))

    for i,main_text_file in enumerate(main_text_files):
        if i%size!=rank:
            continue
        print rank, main_text_file
        tilename = find_tilename(main_text_file)
        epoch_text_file = "{0}/epoch/{1}.epoch.txt".format(indir, tilename)
        out_main = "{0}/main/{1}.fits".format(outdir, tilename)
        out_epoch = "{0}/epoch/{1}.fits".format(outdir, tilename)
        if os.path.exists(out_main) and os.path.exists(out_epoch):
            continue
        try:
            process_text(main_text_file, epoch_text_file, out_main, out_epoch, "r", blind=False, quiet=False, report=report)
        except:
            print "{} did not work".format(out_main)
        if report: 
            return

def mkdir(path):
    try:
        os.makedirs(path)
    except OSError as error:
        if error.errno == errno.EEXIST:
            if os.path.isdir(path):
                pass
            else:
                raise ValueError("Tried to create directory %s but file with name exists already"%path)
        elif error.errno == errno.ENOTDIR: 
            raise ValueError("Tried to create directory %s but some part of the path already exists as a file"%path)
        else:
            raise ValueError("Could not write to %s"%path)

def main():

    description = 'im3shape postprocessor.'
    parser = argparse.ArgumentParser(description=description, add_help=False)
    parser.add_argument('-d', '--indir', type=str, action='store', default="disc")
    parser.add_argument('-o', '--outdir', type=str, action='store', default="fits/disc")

    args = parser.parse_args()

    if not os.path.exists(args.outdir):
        mkdir(args.outdir+"/main")
        mkdir(args.outdir+"/epoch")

    process_dir(args.indir, args.outdir)

main()


