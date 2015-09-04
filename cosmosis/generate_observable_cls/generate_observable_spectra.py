import scipy.interpolate
import numpy as np
import pdb

class Cl_class:
	def load_and_generate_observable_cls(self, block, names, config):

		cl_GG = names.shear_cl_gg
		self.get_l_bins(config)

		# Spectra to use
		shear = config['shear']
		intrinsic_alignments= config['intrinsic_alignments']
		clustering = config['clustering']
		magnification = config['magnification']

		noise = config['noise']
		
		self.Nzbin= block[cl_GG,'nbin']

		# By construction the l sampling should be the same for all the spectra
		self.l = block[cl_GG,'ell']
		Nl = len(self.l)

		self.zbin_edges = [ block['wl_number_density', 'edge_%d'%i] for i in range(1, self.Nzbin + 1) ]

		if shear:
			self.C_ee = np.zeros((self.Nzbin, self.Nzbin, Nl))
			self.C_ee_binned = np.zeros((self.Nzbin, self.Nzbin, self.Nlbin))
		if clustering:
			self.C_nn = np.zeros((self.Nzbin, self.Nzbin, Nl))
			self.C_nn_binned = np.zeros((self.Nzbin, self.Nzbin, self.Nlbin))
		if shear and clustering:
			self.C_ne = np.zeros((self.Nzbin, self.Nzbin, Nl))
			self.C_ne_binned = np.zeros((self.Nzbin, self.Nzbin, self.Nlbin))
	
		for i in range(1, self.Nzbin+1):
			for j in range(1, self.Nzbin+1):
				bin = "bin_%d_%d" %(i,j)
				bin_tr = "bin_%d_%d" %(j,i)

				# The C_GG,II,mm,gg spectra are symmetric
				# This is just bookkeeping to account for the fact we only have half of them
				if (j<i):	a = bin
				else:		a = bin_tr
					
				if shear:
					self.C_ee[i-1][j-1] += block[cl_GG, a]		 						# GG
					if intrinsic_alignments:	
						self.C_ee[i-1][j-1] += block[names.shear_cl_gi, bin] 					# GI
						self.C_ee[i-1][j-1] += block[names.shear_cl_gi, bin_tr]					# IG
						self.C_ee[i-1][j-1] += block[names.shear_cl_ii, a]  					# II

				if clustering:			
					self.C_nn[i-1][j-1] += block['matter_cl', a]							# gg
					if magnification:
						self.C_nn[i-1][j-1] += block[names.galaxy_magnification_cl, bin]			# mg
						self.C_nn[i-1][j-1] += block[names.galaxy_magnification_cl, bin_tr]			# gm
						self.C_nn[i-1][j-1] += block[names.magnification_magnification_cl, a] 			# mm

				if shear and clustering:
					self.C_ne[i-1][j-1] += block[names.ggl_cl, bin]							# gG
					if intrinsic_alignments:
						self.C_ne[i-1][j-1] += block[names.gal_IA_cross_cl, bin]				# gI
						if magnification:
							self.C_ne[i-1][j-1] += block[names.magnification_intrinsic_cl, bin]		# mI
					if magnification:
						self.C_ne[i-1][j-1] += block[names.magnification_shear_cl, bin]				# mG

				if not noise:
					# Finally resample the spectra in the survey angular frequency bins
					if shear:
						self.C_ee_binned[i-1][j-1] = get_binned_cl(self.C_ee[i-1][j-1], self.l, self.lbin_edges )
					if clustering:
						self.C_nn_binned[i-1][j-1] = get_binned_cl(self.C_nn[i-1][j-1], self.l, self.lbin_edges )
					if shear and clustering:
						self.C_ne_binned[i-1][j-1] = get_binned_cl(self.C_ne[i-1][j-1], self.l, self.lbin_edges )

		if noise:
			# Add noise if required	
			self.add_noise(block, shear, clustering)
			# If noise was added earlier, the binning is done here rather than 
			# immediately on loading
			for i in range(1, self.Nzbin+1):
				for j in range(1, self.Nzbin+1):
					if shear:
						self.C_ee_binned[i-1][j-1] = get_binned_cl(self.C_ee[i-1][j-1], self.l, self.lbin_edges )
					if clustering:
						self.C_nn_binned[i-1][j-1] = get_binned_cl(self.C_nn[i-1][j-1], self.l, self.lbin_edges )
					if shear and clustering:
						self.C_ne_binned[i-1][j-1] = get_binned_cl(self.C_ne[i-1][j-1], self.l, self.lbin_edges )

		

	def add_noise(self, block, shear, clustering):
		sur = 'number_density_params'
		sigma_gamma = block[sur, 'shape_dispersion']

		n_binned = get_binned_number_densities(block)

		# Create noise matrices with the same shape as the Cls
		# These are diagonal in the x,z plane (fixed l) and constant along the y axis (constant redshift)
		I = np.identity(block['shear_cl_gg','nbin'])

		N_shape_0 = I * sigma_gamma**2 / (2. * n_binned)
		N_shot_0 = I * 1. / n_binned

		N_shape = [] ; N_shot = []

		for i in range( len(self.C_ee[0][0]) ):
			N_shape += [ N_shape_0 ]
			N_shot += [ N_shot_0 ]

		N_shot = np.swapaxes(N_shot,0,2) ; N_shot = np.swapaxes(N_shot,0,1) 
		N_shape = np.swapaxes(N_shape,0,2) ; N_shape = np.swapaxes(N_shape,0,1)

		# Then add the relevant noise to the Cl matrices
		if shear:	self.C_ee += N_shape
		if clustering:	self.C_nn += N_shot

	def get_l_bins(self, config):
		self.Nlbin = config['Nl']
		self.lmin = config['lmin']
		self.lmax = config['lmax']
		# Define some l bins for this survey
		self.lbin_edges = np.logspace(np.log10(self.lmin), np.log10(self.lmax), self.Nlbin+1)
		self.l_bins = np.exp( (np.log(self.lbin_edges[1:] * self.lbin_edges[:-1]))/2.0 ) 

	def save_cls(self, block, config):

		# Spectra to use
		shear = config['shear']
		intrinsic_alignments= config['intrinsic_alignments']
		clustering = config['clustering']
		magnification = config['magnification']
		
		self.zbin_edges = [ block['wl_number_density', 'edge_%d'%i] for i in range(1, self.Nzbin + 1) ]

		if shear:
			block.put_double_array_1d('galaxy_shape_cl', 'l_bin_edges', self.lbin_edges)
			block.put_double_array_1d('galaxy_shape_cl', 'l_bins', self.l_bins)
			block.put_double_array_1d('galaxy_shape_cl', 'z_bin_edges', self.zbin_edges)
			block.put_int('galaxy_shape_cl', 'nl', self.Nlbin)
			block.put_int('galaxy_shape_cl', 'nz', self.Nzbin)
		if clustering:
			block.put_double_array_1d('galaxy_position_cl', 'l_bin_edges', self.lbin_edges)
			block.put_double_array_1d('galaxy_position_cl', 'l_bins', self.l_bins)
			block.put_double_array_1d('galaxy_position_cl', 'z_bin_edges', self.zbin_edges)
			block.put_int('galaxy_position_cl', 'nl', self.Nlbin)
			block.put_int('galaxy_position_cl', 'nz', self.Nzbin)
		if shear and clustering:
			block.put_double_array_1d('galaxy_position_shape_cross_cl', 'l_bin_edges', self.lbin_edges)
			block.put_double_array_1d('galaxy_position_shape_cross_cl', 'l_bins', self.l_bins)
			block.put_double_array_1d('galaxy_position_shape_cross_cl', 'z_bin_edges', self.zbin_edges)
			block.put_int('galaxy_position_shape_cross_cl', 'nl', self.Nlbin)
			block.put_int('galaxy_position_shape_cross_cl', 'nz', self.Nzbin)
	
		for i in range(1, self.Nzbin+1):
			for j in range(1, self.Nzbin+1):
				bin = "bin_%d_%d" %(i,j)

				if shear:
					block.put_double_array_1d('galaxy_shape_cl', bin, self.C_ee_binned[i-1][j-1])
				if clustering:
					block.put_double_array_1d('galaxy_position_cl', bin, self.C_nn_binned[i-1][j-1])
				if shear and clustering:
					block.put_double_array_1d('galaxy_position_shape_cross_cl', bin, self.C_ne_binned[i-1][j-1])

def get_binned_cl(Cl, l, lbin_edges):
	Cl_binned = np.zeros(len(lbin_edges)-1)
	for i in range(len(lbin_edges)-1):
		lmin = lbin_edges[i] ; lmax = lbin_edges[i+1]
		sel = (l>lmin) & (l<lmax)
		Cl_binned[i] = np.mean(Cl[sel])
	return Cl_binned


def get_binned_number_densities(block):
	num = 'wl_number_density'
	sur = 'number_density_params'

	# Calculate the average number density of galaxies in each redshift bin
	n_binned = []
	nzbin = block[num, 'nbin']					# Number of redshift bins
	neff = ( 180*180 / (np.pi*np.pi) )* block[sur, 'ngal'] 		# Effective number density in galaxies sr^-1
	for i in range(nzbin):
		n_binned += [neff/nzbin]
	n_binned = np.array(n_binned)

	return n_binned


		
