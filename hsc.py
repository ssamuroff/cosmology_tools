import numpy as np
import fitsio as fi
import os


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

	def detect(self, band='r', pointing='(8,7)', config='', weights=True):
		path = '%s/deepCoadd/HSC-%c/9813/%s'%(self.base,band.upper(),pointing)
		full_path = '%s/calexp-HSC-%c-9813-%s.fits.gz'%(path, band.upper(), pointing)
		filename = os.path.basename(full_path)

		print full_path

		if config=='':
			config='%s/seconfig-nersc'%self.base

		template = 'sex %s/deepCoadd/HSC-%c/9813/%s/calexp-HSC-%c-9813-%s.fits.gz -c %s'%(self.base, band.upper(), pointing, band.upper(), pointing, config)

		if weights:
			template += "-WEIGHT_IMAGE %s'[1]'"%filename

		print "running: ", template
		os.system(template)

		print 'Done'



