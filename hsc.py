import numpy as np
import fitsio as fi
import os
import math
import galsim
import copy
import pylab as plt ; plt.switch_backend('agg')

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

	def export_galsim_stamps(self, bands=['r','i','z'], patches=[], mask=False):

		if len(patches)<1:
			patches = patches_all

		igal=0

		if mask:
			suffix="masked"
		else:
			suffix="unmasked"

		for b in bands:

			out_path2 = 'calexp-HSC-%c-9813_galsim_catalogue_%s.fits'%(b.upper(), suffix)
			os.system('rm %s'%out_path2)
			outfile2 = fi.FITS(out_path2, 'rw')

			outdat_all = np.empty(10000000, dtype=[('IDENT', int), ("FLAGS", int), ('RA', float), ('DEC', float), ('GAL_FILENAME', 'S100'), ('GAL_HDU', int)])
			start=0

			for ip,p in enumerate(patches):

				path = '%s/deepCoadd/HSC-%c/9813/%s'%(self.base,b.upper(),p)
				cat_path = '%s/calexp-HSC-%c-9813-%s_cat.fits'%(path, b.upper(), p)
				seg_path = '%s/calexp-HSC-%c-9813-%s_seg.fits'%(path, b.upper(), p)
				coadd_path = '%s/calexp-HSC-%c-9813-%s.fits.gz'%(path, b.upper(), p)
				filename = os.path.basename(cat_path)

				print cat_path

				coadd_data = fi.FITS(coadd_path)['IMAGE'][:,:]
				mask_data = fi.FITS(coadd_path)['MASK'][:,:]
				seg_data = fi.FITS(seg_path)[0][:,:]
				cat_data = fi.FITS(cat_path)[1].read()


				boxsizes = get_boxsizes(cat_data)
				pixel_count = 0

				out_path = '%s/calexp-HSC-%c-9813-%s_galsim_images_%s.fits'%(path, b.upper(), p, suffix)

				outdat = np.zeros(boxsizes.size, dtype=[('IDENT', int), ("FLAGS", int), ('RA', float), ('DEC', float), ('GAL_FILENAME', 'S100'), ('GAL_HDU', int)])
				
				print "Writing cutouts to %s"%out_path
				os.system('rm %s'%out_path)
				outfile = fi.FITS(out_path, 'rw')

				ihdu=0

				for i,row in enumerate(cat_data):
				    x = int(math.floor(row['XWIN_IMAGE']+0.5))
				    y = int(math.floor(row['YWIN_IMAGE']+0.5))

				    boxsize = boxsizes[i]

				    x0 = x-boxsize/2
				    y0 = y-boxsize/2
				    x1 = x+boxsize/2
				    y1 = y+boxsize/2

				    stamp = coadd_data[y0:y1,x0:x1]
				    seg_stamp = seg_data[y0:y1,x0:x1]

				    if not np.sqrt(stamp.size)==boxsize:
				    	im = galsim.ImageD(coadd_data)
				    	bounds = galsim.BoundsI(xmin=x0,xmax=x1,ymin=y0,ymax=y1)

				    	if (stamp.shape[0]==0) or (stamp.shape[1]==0):
				    		continue  


				    	final = np.zeros((boxsize,boxsize))
				    	seg_final = np.zeros((boxsize,boxsize))
				    	dy0 = abs(min(y0-0, 0))
				    	dx0 = abs(min(x0-0, 0))
				    	dy1 = abs(min(coadd_data.shape[0]-y1, 0))
				    	dx1 = abs(min(coadd_data.shape[1]-x1, 0))


				    	final[dy0:boxsize-dy1, dx0:boxsize-dx1]=stamp
				    	seg_final[dy0:boxsize-dy1, dx0:boxsize-dx1]=seg_stamp

				    else:
				        final = stamp
				        seg_final = seg_stamp	


				    if (np.unique(seg_final).size>2):
				    	
				    	if mask:
				    		unmasked_pixels = np.argwhere(seg_final==0)
				    		
				    		sig = np.std(final[:5,])
				    		noise_stamp = np.random.normal(size=final.size).reshape(final.shape) * sig
				    		masked_pixels = (seg_final!=0) & (seg_final!=seg_final[boxsize/2,boxsize/2])
				    		
				    		final[masked_pixels]=noise_stamp[masked_pixels]

				    number = seg_final[boxsize/2,boxsize/2]
				    outfile.write(final)

				    outdat['IDENT'][i] = seg_final[boxsize/2,boxsize/2]
				    outdat['RA'][i] = row['ALPHAWIN_J2000']
				    outdat['DEC'][i] = row['DELTAWIN_J2000']
				    outdat['GAL_FILENAME'][i] = out_path
				    outdat['GAL_HDU'][i] = ihdu
				    outdat['FLAGS'] = cat_data['FLAGS'][cat_data['NUMBER']==number][0]

				    igal+=1
				    ihdu+=1

				#import pdb ; pdb.set_trace()

				for name in outdat_all.dtype.names:
					try:
						outdat_all[name][start:start+outdat[name].size] = outdat[name]
						start+=outdat[name].size
					except:
						import pdb ; pdb.set_trace()

				outfile.close()

			outfile2.write(outdat_all)
			outfile2.close()

		print "Done"


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
				object_data=np.zeros(boxsizes.size, dtype=[('number', int),('start_row', int), ('box_size', int)])
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
				    	im = galsim.ImageD(coadd_data)
				    	bounds = galsim.BoundsI(xmin=x0,xmax=x1,ymin=y0,ymax=y1)

				    	if (stamp.shape[0]==0) or (stamp.shape[1]==0):
				    		continue  


				    	final = np.zeros((boxsize,boxsize))
				    	seg_final = np.zeros((boxsize,boxsize))
				    	dy0 = abs(min(y0-0, 0))
				    	dx0 = abs(min(x0-0, 0))
				    	dy1 = abs(min(coadd_data.shape[0]-y1, 0))
				    	dx1 = abs(min(coadd_data.shape[1]-x1, 0))


				    	final[dy0:boxsize-dy1, dx0:boxsize-dx1]=stamp
				    	seg_final[dy0:boxsize-dy1, dx0:boxsize-dx1]=seg_stamp

				    else:
				        final = stamp
				        seg_final = seg_stamp	

				    image_pixels.append(final.flatten())
				    seg_pixels.append(seg_final.flatten())
				    object_data['number'][i] = seg_final[boxsize/2,boxsize/2]
				    object_data['start_row'][i] = pixel_count
				    object_data['box_size'][i] = boxsize

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


def get_zeropadding(im, bounds):

            # Number of pixels on each axis outside the image bounds
            delta_x_min, delta_x_max, delta_y_min, delta_y_max = 0,0,0,0   
         
            if (bounds.ymin< im.bounds.ymin):
                delta_y_min = im.bounds.ymin - bounds.ymin
                bounds = galsim.BoundsI(ymin=im.bounds.ymin, ymax=bounds.ymax, xmin=bounds.xmin,xmax=bounds.xmax)
            if (bounds.ymax > im.bounds.ymax):
                delta_y_max = bounds.ymax - im.bounds.ymax
                bounds = galsim.BoundsI(ymin=bounds.ymin, ymax=im.bounds.ymax, xmin=bounds.xmin,xmax=bounds.xmax)
            if (bounds.xmin< im.bounds.xmin):
                delta_x_min = im.bounds.xmin - bounds.xmin
                bounds = galsim.BoundsI(ymin=bounds.ymin, ymax=bounds.ymax, xmin=im.bounds.xmin,xmax=bounds.xmax)
            if (bounds.xmax > im.bounds.xmax):
                delta_x_max = bounds.xmax - im.bounds.xmax
                bounds = galsim.BoundsI(ymin=bounds.ymin, ymax=bounds.ymax, xmin=bounds.xmin,xmax=im.bounds.xmax)

            if (bounds.ymin< im.bounds.ymin) or (bounds.ymax > im.bounds.ymax) or (bounds.xmin< im.bounds.xmin) or (bounds.xmax > im.bounds.xmax):
                print 'Stamp bounds are not fully contained by the image bounds'

            return bounds, delta_x_min, delta_x_max, delta_y_min, delta_y_max         


def get_cutout(i, pixels, info):
	i0 = info["start_row"][i]
	b = info["box_size"][i]

	i1 = i0 + b*b

	stamp = pixels[i0:i1].reshape((b,b))

	return stamp





				    










