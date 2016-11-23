import numpy as np
from astropy.table import Table, Column
import tools.shapes as s
import glob, argparse, os

parser = argparse.ArgumentParser(add_help=True)
parser.add_argument('--truth', type=str, default=".", help="Directory in which to look for the truth tables")
parser.add_argument('--output', type=str, default="./hoopoe_y1a1_ra_dec_flags_mag_v2.hdf5", help="Name to save the new object flags under")


args = parser.parse_args()

sim = s.shapecat(truth=args.truth)
sim.load(res=False, truth=True, apply_infocuts=False)

gold_cat_columns = np.dtype([('FLAGS_G', '<i4'), ('FLAGS_R', '<i4'), ('FLAGS_I', '<i4'), ('FLAGS_Z', '<i4'), ('FLAGS_Y', '<i4'), ('RA', '<f8'), ('DEC', '<f8'), ('MAG_AUTO_G', '<f4'), ('MAG_AUTO_R', '<f4'), ('MAG_AUTO_I', '<f4'), ('MAG_AUTO_Z', '<f4'), ('MAG_AUTO_Y', '<f4'), ('DESDM_ZP', '<f8'), ('MODEST', '<f8'), ('HPIX', '<f8'), ('COADD_OBJECTS_ID', '<i8')])

mock_cat = np.zeros(sim.truth.size, dtype=gold_cat_columns)

mock_cat["COADD_OBJECTS_ID"] = sim.truth["DES_id"]
mock_cat["MAG_AUTO_R"] = sim.truth["mag"]
mock_cat["FLAGS_R"] = sim.truth["flags"]
mock_cat["DESDM_ZP"] = sim.truth["cosmos_photoz"]
mock_cat["RA"] = sim.truth["ra"]
mock_cat["DEC"] = sim.truth["dec"]
mock_cat["MODEST"] = sim.truth["star_flag"]+1

mock_cat = Table(mock_cat)

print "Writing mock gold catalogue flags table to %s"%args.output

mock_cat.write(os.path.basename(args.output), path=os.path.dirname(args.output), format="hdf5") 
