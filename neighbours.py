import numpy as np
import astropy.io.fits as pf
import fitsio as fio
import glob
#import treecorr as tc
import meds, fitsio
import galsim
import os, math
import py3shape.i3meds as i3meds
import py3shape.utils as utils
import py3shape.options as i3opt
import tools.diagnostics as di
from scipy.interpolate import Rbf
import tools.arrays as arr
import scipy.spatial as sps
from plots import im3shape_results_plots as i3s_plots
import py3shape as p3s


class neighbours:
    def __init__(self, res=None, truth=None):
    	self.res = res
    	self.truth = truth

    	print "Initialised neighbour catalogue"

    def load_distance_cat(self, catalogue, path="/share/des/disc7/samuroff/des/hoopoe/neighbour_catalogue.fits"):
    	print "Loading neighbour lookup catalogue from %s"%path
    	self.distance_cat = ncat=fi.FITS(path)[1].read(columns=["coadd_objects_id", "neighbour_distance", "truth_index"])

    	dmask =np.in1d(self.distance_cat["coadd_objects_id"],catalogue.res["coadd_objects_id"])
    	self.distance_cat,catalogue.res=di.match_results(self.distance_cat[dmask], catalogue.res)

    def generate_from_truth(self, catalogue, truth_path="/share/des/disc8/cambridge/truth/"):
    	container=np.zeros_like(catalogue.truth)

    	if not hasattr(catalogue, "tiles"):
    		catalogue.tiles = catalogue.res["tilename"].astype("S12")

    	for i,tile in enumerate(catalogue.tiles):
    		print i, tile

    		truth_file = glob.glob("%s/%s*.fz"%(truth_path,tile))
    		if (len(truth_file)==0 ) :
    			continue
    		else:
    			truth_file=truth_file[0]

    		truth = fi.FITS(truth_file)[1].read()
    		select = (hoopoe.res["tilename"].astype("S12")==tile)
    		container[select] = truth[truth["DES_id"]!=-9999][ncat["truth_index"][select]]



    def tune(self, target, ntiles=np.inf, task_number=0, tiles=[]):
    	print "setting up..."

    	self.ids=np.zeros(target.truth.size)
    	self.Rn=np.zeros(target.truth.size)-9999

    	indices_all = np.arange(0, self.truth.size, 1)

    	print "Will process %d tiles"%min(target.tiles.size, ntiles)

    	for i, tile in enumerate(target.tiles):
    		if (i+1)>ntiles:
    			break
    		if len(tiles)>0:
    			if tile not in tiles:
    				continue

    		print i+1, tile
    		select = target.res["tilename"].astype("S12")==tile
    		select_neighbours = (self.truth["DES_id"]<(target.truth["DES_id"][select][-1]+200)) & (self.truth["DES_id"]>(target.truth["DES_id"][select][0]-200))

    		finder = np.arange(0,self.truth["ra"][select_neighbours].size, 1)
    		tile_indices=[]

    		for cid in target.truth["nearest_neighbour_index"][select]:
    			dat = finder[self.truth["DES_id"][select_neighbours]==cid] 
    			if len(dat)>0:
    				tile_indices.append(dat[0])
    			else:
    				tile_indices.append(-1)


    		ids = indices_all[select_neighbours][tile_indices]

    		with open("neighbour_dist-%d.txt"%task_number, "a") as f:
    			out = np.vstack((target.truth["DES_id"][select], ids))
    			np.savetxt(f, out.T)

    	print "Tuned to match source catalogue of %d objects"%len(target.truth)

	def check_tuning(self, target):
		matched = (target==self.matched_to).all()
		if matched:
			print "Yes"
		else:
			print "No"

	def extract_matched_arrays(self, reference, quantity):
		return reference[quantity][self.key]

    def not_so_fast_neighbour_search(self, reference, truth=None, repetition=False):
		if truth is None:
			truth = self.truth_path

		self.tiles = np.unique(reference.res["tilename"]).astype("S12")
		object_tiles = reference.res["tilename"].astype("S12")

		columns=['id', 'DES_id', 'cosmos_ident', 'cosmos_photoz', 'spectral_type', 'flags', 'star_flag', 'nearest_neighbour_pixel_dist', 'nearest_neighbour_index', 'true_g1', 'true_g2', 'intrinsic_e1', 'intrinsic_e2', 'intrinsic_e1_hsm', 'intrinsic_e2_hsm', 'hlr',  'ra', 'dec', 'mag', 'flux', 'nexp', 'mean_psf_e1', 'mean_psf_e2', 'mean_psf_fwhm']

		if not repetition:
			neighbours=[]
		else:
			neighbours=np.zeros_like(self.truth)

		print "Searching for neighbours in %d truth tables."%self.tiles.size

		for i, tile in enumerate(self.tiles):
			
			select = (object_tiles==tile)
			print i+1, tile, object_tiles[select].size

			filename = glob.glob("%s/%s*.fz"%(truth,tile))
			if len(filename)==0:
				print "Truth table is missing."
				continue
			if len(filename)>1:
				print "Warning - multiple truth tables found." 
				print "Will use %s"%filename[0]

			filename = filename[0]

			truth_table = fitsio.FITS(filename)[1].read(columns=columns)
			required = reference.truth["nearest_neighbour_index"][select] 
	
			if repetition:
				finder=np.arange(0,len(truth_table), 1)
				indices=np.array([ finder[truth_table["DES_id"]==ni][0] for ni in required])
				neighbours[select] = truth_table[indices]
			else:
				extract = np.in1d(truth_table["DES_id"],required)  
				neighbours.append(truth_table[extract])

		if repetition:
			neighbours = np.concatenate(neighbours)

		return neighbours



