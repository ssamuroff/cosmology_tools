import numpy as np
import fitsio as fi
import os

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

	def collect_stamps(self, bands=['r','i','z'], patches=[]):

		if len(patches)<1:
			patches = patches_all

		for b in bands:
			for p in patches:

				path = '%s/deepCoadd/HSC-%c/9813/%s'%(self.base,band.upper(),patch)
				cat_path = '%s/calexp-HSC-%c-9813-%s_cat.fits'%(path, band.upper(), patch)
				seg_path = '%s/calexp-HSC-%c-9813-%s_seg.fits'%(path, band.upper(), patch)
				coadd_path = '%s/calexp-HSC-%c-9813-%s.fits.fz'%(path, band.upper(), patch)
				filename = os.path.basename(full_path)

				print full_path

				import pdb ; pdb.set_trace()

				coadd=fi.FITS('')









