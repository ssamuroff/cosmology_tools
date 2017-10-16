import numpy as np
import fitsio as fi
import os
import math

patches_all = ['0,2',  '1,2',  '1,7',  '2,4',  '3,1',  '3,6',  '4,2',  '4,7', '5,3',  '5,8',  '6,4',  '7,0',  '7,5',  '8,2',  '8,7', '0,3',  '1,3',  '2,0',  '2,5',  '3,2',  '3,7',  '4,3',  '4,8',  '5,4',  '6,0',  '6,5',  '7,1',  '7,6',  '8,3', '0,4',  '1,4',  '2,1',  '2,6',  '3,3',  '3,8',  '4,4',  '5,0',  '5,5',  '6,1',  '6,6',  '7,2',  '7,7',  '8,4', '0,5',  '1,5',  '2,2',  '2,7',  '3,4',  '4,0',  '4,5',  '5,1',  '5,6',  '6,2',  '6,7',  '7,3',  '7,8',  '8,5', '1,1',  '1,6',  '2,3',  '3,0',  '3,5',  '4,1',  '4,6',  '5,2',  '5,7',  '6,3',  '6,8',  '7,4',  '8,1',  '8,6']

class dr1:
	def __init__(self, data=''):
		print 'Initialised HSC DR1 data wrapper.'
		if data=='':
			data = '/global/cscratch1/sd/sws/hsc/pdr1_udeep_wide_depth_best'

		self.base = data
		print 'will obtain data from %s'%self.base

		self.data = {}
		return None
	def load_coadd(self, band='r', pointing='(8,7)'):
		path = '%s/deepCoadd/HSC-%c/9813/%s'%(self.base,band.upper(),pointing)
		full_path = '%s/calexp-HSC-%c-9813-%s.fits.gz'%(path, band.upper(), pointing)
		filename = os.path.basename(full_path)

		if filename not in self.data.keys():
			self.data[filename]={}

		self.data[filename]['coadd'] = fi.FITS(filename)
		return 0

	def detect(self, band='r', patch='8,7', config='', weights=True):
		# Work out the paths and make a copy of the coadd
		path = '%s/deepCoadd/HSC-%c/9813/%s'%(self.base,band.upper(),patch)
		full_path = '%s/calexp-HSC-%c-9813-%s.fits.gz'%(path, band.upper(), patch)
		filename = os.path.basename(full_path)

		print full_path

		if config=='':
			config='%s/seconfig-nersc'%self.base

		os.system('cp %s tmp.fits.gz'%full_path)
		os.system('gunzip tmp.fits.gz')

		# Construct the command
		template = "sex tmp.fits'[1]' -c %s"%(config)

		if weights:
			template += " -WEIGHT_IMAGE tmp.fits'[3]'"

		# Run the command
		print "running: ", template
		os.system(template)

		# Now clean up
		os.system('rm tmp.fits.gz')
		os.system('rm tmp.fits')

		os.system('mv cat.fits %s/calexp-HSC-%c-9813-%s_cat.fits'%(path, band.upper(), patch))
		os.system('mv seg.fits %s/calexp-HSC-%c-9813-%s_seg.fits'%(path, band.upper(), patch))

		print 'Done'

	def bulk_detect(self, bands=['r','i','z'], patches=[]):
		if len(patches)<1:
			patches = patches_all

		for b in bands:
			for p in patches:
				self.detect(band=b, patch=p)

		print "Done all pointings requested"
		return None


	def collect_stamps(self, bands=['r','i','z'], patches=[], mask=False):

		if len(patches)<1:
			patches = patches_all

		for b in bands:
			for p in patches:

				path = '%s/deepCoadd/HSC-%c/9813/%s'%(self.base,b.upper(),p)
				cat_path = '%s/calexp-HSC-%c-9813-%s_cat.fits'%(path, b.upper(), p)
				seg_path = '%s/calexp-HSC-%c-9813-%s_seg.fits'%(path, b.upper(), p)
				coadd_path = '%s/calexp-HSC-%c-9813-%s.fits.gz'%(path, b.upper(), p)
				filename = os.path.basename(cat_path)

				print cat_path

				coadd_data = fi.FITS(coadd_path)['IMAGE'][:,:]
				seg_data = fi.FITS(seg_path)[0][:,:]
				cat_data = fi.FITS(cat_path)[1].read()


				boxsizes = get_boxsizes(cat_data)

				image_pixels = []
				seg_pixels = []
				object_data=np.zeros(boxsizes.size, dtype=[('number', int),('start_row', int)])
				pixel_count = 0

				for i,row in enumerate(cat_data):
				    x = int(math.floor(row['XWIN_IMAGE']+0.5))
				    y = int(math.floor(row['YWIN_IMAGE']+0.5))

				    boxsize = boxsizes[i]
				    pixel_count+=boxsize*boxsize

				    x0 = x-boxsize/2
				    y0 = y-boxsize/2
				    x1 = x+boxsize/2
				    y1 = y+boxsize/2

				    stamp = coadd_data[y0:y1,x0:x1]
				    seg_stamp = seg_data[y0:y1,x0:x1]

				    if not np.sqrt(stamp.size)==boxsize:
				    	continue

				    image_pixels.append(stamp.flatten())
				    seg_pixels.append(seg_stamp.flatten())
				    object_data['number'][i] = seg_stamp[boxsize/2,boxsize/2]
				    object_data['start_row'][i] = pixel_count

				    if mask & np.unique(seg_stamp).size>2:
				    	import pdb ; pdb.set_trace()

				meds_path = '%s/calexp-HSC-%c-9813-%s_meds.fits'%(path, b.upper(), p)
				print "Writing cutouts to %s"%meds_path
				os.system('rm %s'%meds_path)
				meds = fi.FITS(meds_path, 'rw')

				meds.write(np.concatenate(image_pixels))
				meds[-1].write_key('EXTNAME','image_cutouts')
				meds.write(np.concatenate(seg_pixels))
				meds[-1].write_key('EXTNAME','seg_cutouts')
				meds.write(object_data)
				meds[-1].write_key('EXTNAME','object_data')
				meds.close()

		print "Done"






def get_boxsizes(cat_data):
    sig=2*cat_data['FLUX_RADIUS']/2.3548
    epsilon= 1- cat_data['B_WORLD']/cat_data['A_WORLD']
    boxsize = 2*5*sig*(1+epsilon)

    min_boxsize, max_boxsize = 32, 256

    isarray = True
    boxsize = (boxsize+0.5).astype(int)

    # Replace any entries with infinite or overly large or small boxes
    np.putmask(boxsize, np.invert(np.isfinite(boxsize)), 32)
    np.putmask(boxsize, boxsize<min_boxsize, min_boxsize)
    np.putmask(boxsize, boxsize>max_boxsize, max_boxsize)

    # Round up to a factor of 2**N or 3*2**N
    allowed_boxsizes = 2**np.arange(0,1 + np.log(max_boxsize)/np.log(2),1 )
    allowed_boxsizes = np.concatenate((allowed_boxsizes, 3*allowed_boxsizes))
    allowed_boxsizes.sort()
    allowed_boxsizes=allowed_boxsizes[(allowed_boxsizes>=min_boxsize) & (allowed_boxsizes<=max_boxsize) ]

    for i, bs in enumerate(boxsize):
        if bs not in allowed_boxsizes:
            boxsize[i] = allowed_boxsizes[allowed_boxsizes>bs][0] 
        else:
            continue

    return boxsize  




				    










