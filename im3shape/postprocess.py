"""

This is my attmept to combine all the operations using astropy to make everything simpler.

"""
import glob
import os
#import i3astro
import sys
import argparse
from . import cut_data
from . import blind_catalog
from . import split_cat
from . import append_wcs
from . import extra_columns
from . import db
from . import merge_bulge_disc
import astropy.table
import numpy as np
import re
import errno

class FakeCommunicator(object):
    def Get_rank(self):
        return 0
    def Get_size(self):
        return 1
    def Barrier(self):
        pass

def read_split_table(filenames, fatal_errors):
    tables = []
    for filename in filenames:
        try:
            t = astropy.table.Table.read(filename, format='ascii')
            newdt=[]
            for n in t.dtype.names:
                if ("covmat" in n) or ("hsm" in n) or  ("model" in n) or ("min" in n) or ("max" in n) or ("levmar" in n):
                     newdt.append('float64')
                else:
                     newdt.append( t.dtype[n].name.replace("string", "S") )

            tables.append(astropy.table.Table(np.array(t),names=t.dtype.names,dtype=tuple(newdt) ))
             
        except:
            if fatal_errors:
                raise
            else:
                print "FILE BROKEN/EMPTY: ", filename
    if not tables:
        return None
    table = astropy.table.vstack(tables, join_type='exact', metadata_conflicts='silent')
    return table

def postprocess_hierarchical(job_name, band, comm, blind=True, quiet=False, fatal_errors=False):
    #files will be in {job_name}/{model}/{meds_name}/{split_num}.{main_epoch}.txt
    #load all the main files and epoch files
    rank = comm.Get_rank()
    size = comm.Get_size()
    if rank==0:
        if blind:
            print "BEGINNING POSTPROCESS"
        else:
            print "BEGINNING POSTPROCESS (BLINDING IS OFF!!)"

    disc_indir = "{0}/disc".format(job_name)
    bulge_indir = "{0}/bulge".format(job_name)
    disc_dir = "{0}/disc-fits".format(job_name)
    bulge_dir = "{0}/bulge-fits".format(job_name)
    bord_dir = "{0}/bord-fits".format(job_name)
    disc_main_dir = "{0}/main".format(disc_dir)
    bulge_main_dir = "{0}/main".format(bulge_dir)
    bord_main_dir = "{0}/main".format(bord_dir)
    disc_epoch_dir = "{0}/epoch".format(disc_dir)
    bulge_epoch_dir = "{0}/epoch".format(bulge_dir)
    bord_epoch_dir = "{0}/epoch".format(bord_dir)

    
    dirs = [disc_dir, bulge_dir, bord_dir, 
            disc_main_dir, bulge_main_dir, bord_main_dir, 
            disc_epoch_dir, bulge_epoch_dir, bord_epoch_dir]

    for dirname in dirs:
        mkdir(dirname)

    for in_dir, main_dir, epoch_dir in zip(
            [bulge_indir, disc_indir], 
            [bulge_main_dir, disc_main_dir], 
            [bulge_epoch_dir, disc_epoch_dir]):

        dir_names = glob.glob("{0}/DES*".format(in_dir))
        dir_names.sort()
        dir_names_incomplete = []
        existing_tiles = set()
        
        for dir_name in dir_names:
            tilename = find_tilename(dir_name)
            main_name = "{0}/{1}.fits".format(main_dir, tilename)
            epoch_name = "{0}/{1}.fits".format(epoch_dir, tilename)
            if os.path.exists(main_name) and os.path.exists(epoch_name):
                if rank==0:
                    print "Skipping completed tile", tilename
            elif tilename in existing_tiles:
                if rank==0:
                    print "Skipping repeated tile", tilename
            else:
                dir_names_incomplete.append(dir_name)
            existing_tiles.add(tilename)

        print "Found {} tiles left to process for {}".format(len(dir_names_incomplete), in_dir)

        for i,dir_name in enumerate(dir_names_incomplete):
            if i%size!=rank: continue
            tilename = find_tilename(dir_name)
            main_name = "{0}/{1}.fits".format(main_dir, tilename)
            epoch_name = "{0}/{1}.fits".format(epoch_dir, tilename)

            print "Processor {} working on {} from {}".format(rank, tilename, in_dir)

            #Find all the split tables from the various jobs
            main_files = glob.glob("{0}/*.main.txt".format(dir_name))
            epoch_files = glob.glob("{0}/*.epoch.txt".format(dir_name))
            print "{0}/{1}/*.main.txt".format(in_dir, dir_name)
            if not quiet:
                print "Found {} main files and {} epoch files for {}".format(len(main_files), len(epoch_files), tilename)
            
            #Read them all in
            main_table = read_split_table(main_files, fatal_errors)
            epoch_table = read_split_table(epoch_files, fatal_errors)

            if main_table is None:
                print "No valid main files for {} {}".format(in_dir, tilename)
                continue
            if epoch_table is None:
                print "No valid epoch files for {} {}".format(in_dir, tilename)
                continue

            #Postprocess with the usual stuff
            try:
                main_table, epoch_table = process_table(main_table, epoch_table, band, blind=blind, quiet=quiet)
            except:
                if fatal_errors:
                    raise
                else:
                    print "FAILED TO POSTPROCESS", tilename

            #
            main_table.write(main_name)
            epoch_table.write(epoch_name)

    comm.Barrier()
    if rank==0:
        print "BEGINNING MERGE!"

    #Need to find all the tiles where we have all four of the files
    #that we need
    bulge_main_files = glob.glob("{0}/*.fits".format(bulge_main_dir))
    disc_main_files = glob.glob("{0}/*.fits".format(disc_main_dir))
    bulge_epoch_files = glob.glob("{0}/*.fits".format(bulge_epoch_dir))
    disc_epoch_files = glob.glob("{0}/*.fits".format(disc_epoch_dir))

    #Find the tiles which have files in all the different directories
    completed_tiles  = {find_tilename(f) for f in bulge_main_files}
    completed_tiles.intersection_update({find_tilename(f) for f in  disc_main_files})
    completed_tiles.intersection_update({find_tilename(f) for f in  bulge_epoch_files})
    completed_tiles.intersection_update({find_tilename(f) for f in  disc_epoch_files})

    for t in [bulge_main_files, disc_main_files, bulge_epoch_files, disc_epoch_files]:
        t.sort()

    bulge_main_files = [f for f in bulge_main_files if find_tilename(f) in completed_tiles]
    disc_main_files = [f for f in disc_main_files if find_tilename(f) in completed_tiles]
    bulge_epoch_files = [f for f in bulge_epoch_files if find_tilename(f) in completed_tiles]
    disc_epoch_files = [f for f in disc_epoch_files if find_tilename(f) in completed_tiles]

    print "Found complete information for {} tiles".format(len(completed_tiles))

    for i, (bulge_main, bulge_epoch, disc_main, disc_epoch) in enumerate(zip(
        bulge_main_files, bulge_epoch_files, disc_main_files, disc_epoch_files)):

        if i%size!=rank: continue

        tilename = find_tilename(bulge_main)
        print "MERGING ", tilename
        bord_main = "{0}/{1}.fits".format(bord_main_dir, tilename)
        bord_epoch = "{0}/{1}.fits".format(bord_epoch_dir, tilename)

        try:
            merge_bulge_disc.merge(bulge_main, bulge_epoch, disc_main, disc_epoch, bord_main, bord_epoch)
        except:
            if fatal_errors:
                raise
            else:
                print "FAILED TO MERGE", tilename


    #Now merge bulge and disc


def process_multi_tile_text(main_file, epoch_file, outdir, band, blind=True, quiet=True):
    mkdir(outdir+"/text")
    mkdir(outdir+"/text/main")
    mkdir(outdir+"/text/epoch")
    mkdir(outdir+"/fits/")
    mkdir(outdir+"/fits/main")
    mkdir(outdir+"/fits/epoch")


    filenames = split_by_tile(main_file, epoch_file, outdir+"/text")

    for n,(tilename, main_file, epoch_file) in enumerate(filenames):
        print "{0} - Working on {0}".format(n,tilename)
        table = astropy.table.Table.read(main_file, format='ascii')
        epochs = astropy.table.Table.read(epoch_file, format='ascii')
        table, epochs = process_table(table, epochs, band, blind=blind,quiet=quiet)     
        out_main = outdir+"/fits/main/"+tilename+".txt"
        out_epoch = outdir+"/fits/epoch/"+tilename+".txt"
        print 'Saving to FITS'
        table.write(out_main, format='fits')
        epochs.write(out_epoch, format='fits')

def process_text(main_file, epoch_file, out_main, out_epoch, band, blind=True,quiet=True, report=False):
    if os.path.exists(out_main):
        print 'Skipping sub-file', out_main
        return
    table = astropy.table.Table.read(main_file, format='ascii')
    epochs = astropy.table.Table.read(epoch_file, format='ascii')

    if len(table)==0:
        print "No rows found for ", outfile
        return

    table, epochs = process_table(table, epochs, band, blind=blind,quiet=quiet,report=report)
    if report:
        print "Not saving output file as we are just reporting cuts"
        return

    print 'Saving to FITS'
    table.write(out_main, format='fits')
    epochs.write(out_epoch, format='fits')

def process_db(connection, name, tile, outfile, epochfile, band, blind=True, quiet=True):

    sql="select * from {0}_main where tilename='{1}'".format(name, tile)
    print sql
    table = db.table_from_sql(connection, sql)
    if table is None:
        print "NO OBJECTS FOR ", name, tile
        return
    print "Loaded main table"
    sql="select * from {0}_epoch where tilename='{1}'".format(name, tile)
    print sql
    epochs = db.table_from_sql(connection, sql)
    print "Loaded epoch table"
    table, epochs = process_table(table, epochs, band, blind=blind,quiet=quiet)

    print 'Saving to FITS'
    table.write(outfile, format='fits')
    epochs.write(epochfile, format='fits')


def lowercase(table):
    for col in table.colnames[:]:
        if col!=col.lower():
            table.rename_column(col, col.lower())


def process_table(table, epochs, band, blind=True,quiet=True, report=False):
    print 'Removing large boxes'
    too_large_rows=np.where(table['stamp_size']>100)[0]
    table.remove_rows(too_large_rows)
    
    print "Lower casing"
    lowercase(table)
    lowercase(epochs)

    print "Sorting and grouping"
    epochs.sort(["id", "expnum"])
    table.sort("identifier")

    print 'Adding chi2_pixel'
    chi2 = -2 * table['likelihood']
    effective_npix = table['stamp_size']**2 * table['n_exposure'] * (1-table['mean_mask_fraction'])
    chi2_pixel = chi2 / effective_npix
    col = astropy.table.Column(name='chi2_pixel', data=chi2_pixel)
    table.add_column(col)

    print 'Renaming identifier->coadd_objects_ids'
    table.rename_column("identifier", "coadd_objects_id")
    epochs.rename_column('id', 'coadd_objects_id')

    #print 'NOT NOT NOT NOT NOT NOT NOT Adding extra columns flags'
    print "adding new columns"
    extra_columns.add_standard_cols(table)

    print 'Flagging'
    #error cuts
    cuts = cut_data.error_cuts + getattr(cut_data, "error_cuts_%s"%band)
    flags, fractions = cut_data.compute_flags(table, cuts, verb=not quiet)
    col = astropy.table.Column(name='error_flag', data=flags)
    table.add_column(col)
    if report:
        print "INFO CUTS"
        cut_data.report(cuts)


    #info cuts
    cuts = cut_data.info_cuts1 + getattr(cut_data, "info_cuts_%s"%band) + cut_data.info_cuts2
    flags, fractions = cut_data.compute_flags(table, cuts, verb=not quiet)
    col = astropy.table.Column(name='info_flag', data=flags)
    table.add_column(col)

    if report:
        print "ERROR CUTS"
        cut_data.report(cuts)
        return table, epochs

    #blinding
    if blind:
        print 'Blinding'
        e1,e2 = blind_catalog.blind_arrays(table['e1'],table['e2'])
        table['e1'][:] = e1
        table['e2'][:] = e2

    print "Reticulating splines"
    return table, epochs


tile_pattern = re.compile(r'DES[0-9][0-9][0-9][0-9][+-][0-9][0-9][0-9][0-9]')

def find_tilename(name):
    m = tile_pattern.search(name)
    if m is None:
        return "unknown"
    return m.group()


def find_file_pairs(main_dir, epoch_dir):
    main_files = glob.glob(main_dir+"/*.txt")
    epoch_files = glob.glob(epoch_dir+"/*.txt")
    main_tiles = {}
    epoch_tiles = {}
    for filename in main_files:
        tile = find_tilename(filename)
        main_tiles[tile] = filename

    output = []
    for epoch_file in epoch_files:
        tile = find_tilename(epoch_file)
        main_file = main_tiles[tile]
        output.append((main_file, epoch_file, tile))
    return output

def mkdir(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

def split_by_tile(main_file, epoch_file, outdir):
    #Loop through the main file, copying to the correct
    f = open(main_file)
    header_line = f.next()
    cols = header_line.strip('#').split()
    tilename_index = cols.index("tilename")
    identifier_index = cols.index("identifier")
    outfiles = {}
    tilenames = {}
    print "Splitting multi-tile file {0} into dir {1}".format(main_file, outdir)
    for line in f:
        words = line.split()
        tilename = words[tilename_index]
        identifier = int(words[identifier_index])
        tile_file = outfiles.get(tilename)
        if tile_file is None:
            tile_file = open(outdir + "/main/"+tilename+".txt", "w")
            outfiles[tilename] = tile_file
            tile_file.write(header_line)
        tile_file.write(line)
        tilenames[identifier] = tilename
    f.close()
    for tile_file in outfiles.values():
        tile_file.close()
    #Now the epoch file - the tilename is not in there so
    #we look it up in the dict we just made
    f = open(epoch_file)
    header_line = f.next()
    cols = header_line.strip('#').split()
    identifier_index = cols.index("ID")
    epoch_outfiles = {}

    print "Found {0} objects in multi-tile file".format(len(tilenames))
    print "Splitting multi-tile epoch file {0} into dir {1}".format(epoch_file, outdir)

    for line in f:
        words = line.split()
        identifier = int(words[identifier_index])
        tilename = tilenames[identifier]
        tile_file = epoch_outfiles.get(tilename)
        if tile_file is None:
            tile_file = open(outdir + "/epoch/"+tilename+".txt", "w")
            epoch_outfiles[tilename] = tile_file
            tile_file.write(header_line)
        tile_file.write(line)
    f.close()

    for tile_file in epoch_outfiles.values():
        tile_file.close()


    tilenames = outfiles.keys()
    filenames = [(tilename, outfiles[tilename].name, epoch_outfiles[tilename].name) for tilename in tilenames]
    return filenames


dt=np.dtype([('ra_as', '<f8'), ('dec_as', '<f8'), ('e1', '<f8'), ('e2', '<f8'), ('radius', '<f8'), ('radius_ratio', '<f8'), ('bulge_A', '<f8'), ('disc_A', '<f8'), ('bulge_index', '<f8'), ('disc_index', '<f8'), ('identifier', '<i8'), ('time', '<f8'), ('bulge_flux', '<f8'), ('disc_flux', '<f8'), ('flux_ratio', '<f8'), ('snr', '<f8'), ('old_snr', '<f8'), ('min_residuals', '<f8'), ('max_residuals', '<f8'), ('model_min', '<f8'), ('model_max', '<f8'), ('likelihood', '<f8'), ('levmar_start_error', '<f8'), ('levmar_end_error', '<f8'), ('levmar_resid_grad', '<f8'), ('levmar_vector_diff', '<f8'), ('levmar_error_diff', '<f8'), ('levmar_comp_grad', '<f8'), ('levmar_iterations', '<i8'), ('levmar_reason', '<i8'), ('levmar_like_evals', '<i8'), ('levmar_grad_evals', '<i8'), ('levmar_sys_evals', '<i8'), ('mean_flux', '<f8'), ('number_varied_params', '<i8'), ('covmat_0_0', '<f8'), ('covmat_0_1', '<f8'), ('covmat_0_2', '<f8'), ('covmat_0_3', '<f8'), ('covmat_0_4', '<f8'), ('covmat_0_5', '<f8'), ('covmat_1_0', '<f8'), ('covmat_1_1', '<f8'), ('covmat_1_2', '<f8'), ('covmat_1_3', '<f8'), ('covmat_1_4', '<f8'), ('covmat_1_5', '<f8'), ('covmat_2_0', '<f8'), ('covmat_2_1', '<f8'), ('covmat_2_2', '<f8'), ('covmat_2_3', '<f8'), ('covmat_2_4', '<f8'), ('covmat_2_5', '<f8'), ('covmat_3_0', '<f8'), ('covmat_3_1', '<f8'), ('covmat_3_2', '<f8'), ('covmat_3_3', '<f8'), ('covmat_3_4', '<f8'), ('covmat_3_5', '<f8'), ('covmat_4_0', '<f8'), ('covmat_4_1', '<f8'), ('covmat_4_2', '<f8'), ('covmat_4_3', '<f8'), ('covmat_4_4', '<f8'), ('covmat_4_5', '<f8'), ('covmat_5_0', '<f8'), ('covmat_5_1', '<f8'), ('covmat_5_2', '<f8'), ('covmat_5_3', '<f8'), ('covmat_5_4', '<f8'), ('covmat_5_5', '<f8'), ('expnum', '<i8'), ('nparam_varied', '<i8'), ('stamp_size', '<i8'), ('mean_rgpp_rp', '<f8'), ('fails_rgpp_rp', '<i8'), ('mean_psf_e1_sky', '<f8'), ('fails_psf_e1_sky', '<i8'), ('mean_psf_e2_sky', '<f8'), ('fails_psf_e2_sky', '<i8'), ('mean_psf_fwhm', '<f8'), ('fails_psf_fwhm', '<i8'), ('mean_unmasked_flux_frac', '<f8'), ('fails_unmasked_flux_frac', '<i8'), ('mean_model_edge_mu', '<f8'), ('fails_model_edge_mu', '<i8'), ('mean_model_edge_sigma', '<f8'), ('fails_model_edge_sigma', '<i8'), ('mean_edge_mu', '<f8'), ('fails_edge_mu', '<i8'), ('mean_edge_sigma', '<f8'), ('fails_edge_sigma', '<i8'), ('mean_hsm_psf_e1_sky', '<f8'), ('fails_hsm_psf_e1_sky', '<i8'), ('mean_hsm_psf_e2_sky', '<f8'), ('fails_hsm_psf_e2_sky', '<i8'), ('mean_hsm_psf_sigma', '<f8'), ('fails_hsm_psf_sigma', '<i8'), ('mean_hsm_psf_rho4', '<f8'), ('fails_hsm_psf_rho4', '<i8'), ('mean_mask_fraction', '<f8'), ('fails_mask_fraction', '<i8'), ('round_snr', '<f8'), ('fails_round_snr', '<i8'), ('round_snr_mw', '<f8'), ('fails_round_snr_mw', '<i8'), ('bands', 'S5'), ('dec', '<f8'), ('tilename', 'S12'), ('ra', '<f8')])


parser = argparse.ArgumentParser(description="Turn an ascii catalog into fits")
parser.add_argument('-k', '--skip', action='store_true', default=False, help="Skip tiles that are already done and have FITS files there") 
parser.add_argument('-v', '--verbose', action='store_true', default=False, help="Print out each cat file as loaded")
parser.add_argument('--no-blinding', default=False, action='store_true',help="Whether to blind the results")
parser.add_argument('--fatal-errors', default=False, action='store_true',help="Errors reading any files are fatal. For debugging.")
parser.add_argument('--mpi', default=False, action='store_true',help="Split by MPI process.")
parser.add_argument("job_name", help="Base directory for run")
parser.add_argument("band", help="The band. Obvs.")


if __name__=="__main__":
    args = parser.parse_args()
    if args.mpi:
        import mpi4py.MPI
        comm = mpi4py.MPI.COMM_WORLD
    else:
        comm = FakeCommunicator()

    postprocess_hierarchical(args.job_name, args.band, comm, blind=not args.no_blinding, quiet=not args.verbose, fatal_errors=args.fatal_errors)

