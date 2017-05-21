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
import scipy.interpolate as spi
import tools.arrays as arr
import tools.neighbours as nf
from plots import im3shape_results_plots as i3s_plots
import py3shape as p3s

class calculator:
	def __init__(self, image, weight="tophat", wparams=[]):
		self.nx,self.ny = image.shape
		print "Initialising moments calculator with %dx%d image"%image.shape
		self.image = image
		self.x = np.linspace(0,self.nx-1,self.nx) - self.nx/2
		self.y = np.linspace(0,self.ny-1,self.ny) - self.ny/2
		self.xx,self.yy = np.meshgrid(self.x,self.y)

		self.setup_weight(wtype=weight, params=wparams)

	def setup_weight(self, wtype="tophat", params=[], normalise=True):
		if (wtype=="tophat"):
			wgrid = np.ones_like(self.image)
		elif (wtype=="gaussian"):
			# mu is a single number
			# cov is a 2x2 covariance matrix
			mu, cov = params
			wgrid = np.exp(-( (self.xx-mu)*(self.xx-mu)/2*cov[0,0] + (self.yy-mu)*(self.yy-mu)/2*cov[1,1] - (self.yy-mu)*(self.xx-mu)/2*cov[0,1] ) )

		if normalise:
			norm = np.trapz(wgrid, self.x, axis=1)
			norm = np.trapz(norm, self.y)
			self.wt = wgrid/norm
		else:
			self.wt = wgrid

	def upsample(self, factor=5):

		interpolator = spi.interp2d(self.x, self.y, self.image)
		self.xf, self.yf = np.linspace(0,self.nx-1,self.nx*factor) - self.nx/2, np.linspace(0,self.ny-1,self.ny*factor) - self.ny/2
		self.uimage =  interpolator(self.xf, self.yf)
		self.xxf, self.yyf = np.meshgrid(self.xf,self.yf)

	def find_centroid(self):
		norm = np.trapz(self.image*self.wt, self.x, axis=1)
		norm = np.trapz(norm, self.y)

		x0 = np.trapz(self.image*self.xx*self.wt, self.x, axis=1)
		x0 = np.trapz(x0, self.y)/norm

		y0 = np.trapz(self.image*self.yy*self.wt, self.x, axis=1)
		y0 = np.trapz(y0, self.y)/norm

		print x0,y0

		self.centroid = {"x":x0,"y":y0}
		self.norm = norm

	def find_moment(self, mtype=("x","x")):
		kernels = []

		for axis in mtype:
			coord = getattr(self,axis*2)
			c0 = self.centroid[axis]
			kernels.append(coord-c0)

		integrand = self.image * kernels[0] * kernels[1]

		q = np.trapz(integrand, self.x, axis=1)
		q = np.trapz(q, self.y, axis=0)

		return q

	def find_moments(self):
		xx = self.find_moment(mtype=("x","x"))
		yy = self.find_moment(mtype=("y","y"))
		xy = self.find_moment(mtype=("x","y"))
		return [xx,yy,xy]






