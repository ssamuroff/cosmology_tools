import scipy.interpolate
import numpy as np
import pdb

class Cl_class:
	def __init__(self, block, config):
		# Decide which spectra to use
		shear = config['shear']
		clustering = config['clustering']

		self.types = []
		if shear:	
			self.types += [('ee', 'galaxy_shape_cl')]
			self.C_ee={}
		if clustering:	
			self.types += [('nn', 'galaxy_position_cl')]
			self.C_nn={}
		if shear and clustering:
			self.types += [('ne', 'galaxy_position_shape_cross_cl')] 
			self.C_ne={}

		# Survey parameters
		self.lbin_edges = block[self.types[0][1],'l_bin_edges']
		self.l_bins = block[self.types[0][1],'l_bins']
		self.dl = self.lbin_edges[1:] - self.lbin_edges[:-1]
		
		self.Nlbin = len(self.lbin_edges)-1 

		self.Nzbin= block['wl_number_density','nbin']
		self.zbin_edges = [ block['wl_number_density', 'edge_%d'%i] for i in range(1, self.Nzbin + 1) ]

		self.A = ( np.pi*np.pi / (180*180) ) * block['number_density_params', 'survey_area'] # Survey area in sr

	def load_Cls(self, block):

		for i in range(1,self.Nzbin+1):
			for j in range(1,self.Nzbin+1):
				bin = "bin_%d_%d" %(i,j)
				
				for spect in self.types:
					typ = spect[0]
					section = spect[1]
					load_string = 'self.C_%s[bin] = block[section, bin]' %typ 
					exec load_string 

	def calculate_all(self, block):
		print 'Evaluating covariance matrices'
		
		for spect1 in self.types:
			for spect2 in self.types:
				self.compute_covariance_matrix(block, mode=[spect1[0],spect2[0]])
		block.put_double_array_1d('shear_covariance', 'l_bin_edges', self.lbin_edges)
		block.put_double_array_1d('shear_covariance', 'z_bin_edges', self.zbin_edges)

	def compute_covariance_matrix(self, block, mode):

		# The covariance involves various cross terms so the four 
		# indices need to be extracted, two from each spectrum
		indices = [('alpha_gamma', mode[0][0], mode[1][1]), 
			   ('beta_delta', mode[0][1], mode[1][0]), 
			   ('alpha_delta', mode[0][0], mode[1][0]), 
			   ('beta_gamma', mode[0][1], mode[1][1])]

		cov = np.zeros((self.Nzbin, self.Nzbin, self.Nzbin, self.Nzbin, self.Nlbin))

		p = 2. * np.pi / ( (2 * self.l_bins + 1) * self.dl * self.A ) 

		# Set up the four spectra needed for this covariance matrix
		# The point of the following is that the en spectra are not
		# saved separately, since C^ij_ne = C^ji_en 
		r_alpha_gamma = False ; r_beta_delta = False ; r_alpha_delta = False ; r_beta_gamma = False

		for sp in indices:
			try:
				exec 'C_%s = self.C_%c%c ; ' %(sp[0], sp[1], sp[2])
			except:
				exec 'C_%s = self.C_%c%c ; r_%s = True' %(sp[0], sp[2], sp[1], sp[0])

		#pdb.set_trace()
		
		# Then loop over all combinations of the four redshift bins
		for i in range(1,self.Nzbin+1):
			for j in range(1,self.Nzbin+1):
				for k in range(1,self.Nzbin+1):
					for m in range(1,self.Nzbin+1):
						a0= i ; b0 = k  ;  a1= j ; b1 = m
						a2= i ; b2 = m  ;  a3= j ; b3 = k
						if r_alpha_gamma:	a0,b0 = b0,a0
						if r_beta_delta:	a1,b1 = b1,a1
						if r_alpha_delta:	a2,b2 = b2,a2
						if r_beta_gamma:	a3,b3 = b3,a3
						Cl0 = C_alpha_gamma['bin_%d_%d' %(a0,b0)]
						Cl1 = C_beta_delta['bin_%d_%d' %(a1,b1)]
						Cl2 = C_alpha_delta['bin_%d_%d' %(a2,b2)]
						Cl3 = C_beta_gamma['bin_%d_%d' %(a3,b3)]
						cov[i-1][j-1][k-1][m-1] = p * (Cl0 * Cl1 + Cl2 * Cl3)

		save_cov = 'self.cov_%s%s = cov' %(mode[0],mode[1])
		exec save_cov

		block.put_double_array_nd('shear_covariance', 'Cov_%s%s'%(mode[0],mode[1]), cov)
		
		print 'Saved covariance matrix %s%s' %(mode[0],mode[1])

	def plot_cov(self,mode):
		import matplotlib.pyplot as plt

		exec 'cov= self.cov_%c%c%c%c'%(mode[0],mode[1],mode[2],mode[3])
		plt.subplot(2,2,1)
		plt.imshow(np.diag(cov[0][0][0][0]), interpolation='None', origin='lower')
		plt.colorbar() 

		plt.subplot(2,2,2)
		plt.imshow(np.diag(cov[0][1][0][0]), interpolation='None', origin='lower')
		plt.colorbar() 

		plt.subplot(2,2,3)
		plt.imshow(np.diag(cov[0][0][0][1]), interpolation='None', origin='lower')
		plt.colorbar() 

		plt.subplot(2,2,4)
		plt.imshow(np.diag(cov[0][1][0][1]), interpolation='None', origin='lower')
		plt.colorbar()
		plt.show() 


	def output_binned_cls(self, Nl, lmin, lmax):
		# Define some l bins, evenly separated in log space
		l = self.C_GG['l']
		lbin_edges = np.logspace(np.log10(lmin), np.log10(lmax), Nl+1)

		dl = lbin_edges[1:]-lbin_edges[:-1]

		l_binned = (lbin_edges[1:]+lbin_edges[:-1])/2.0
		
		Cl_ee_binned = np.zeros(( self.Nzbin, self.Nzbin, Nl ))
		Cl_nn_binned = np.zeros(( self.Nzbin, self.Nzbin, Nl ))
		Cl_ne_binned = np.zeros(( self.Nzbin, self.Nzbin, Nl ))
		# Then loop over all combinations of the two redshift bins
		for i in range(self.Nzbin):
			for j in range(self.Nzbin):
				# Resample the theory spectra in the l bins defined
				Cl_ee_binned[i][j] = get_binned_cl(self.C_ee[i][j], l, lbin_edges)
				Cl_nn_binned[i][j] = get_binned_cl(self.C_nn[i][j], l, lbin_edges)
				Cl_ne_binned[i][j] = get_binned_cl(self.C_ne[i][j], l, lbin_edges)

		out_data = {'Cl_ee': Cl_ee_binned, 'Cl_nn': Cl_nn_binned, 'Cl_ne': Cl_ne_binned, 'l_bin_edges': lbin_edges, 'z_bin_edges': self.zbin_edges}

		import cPickle as pickle
		fil = open('cl_datavector.p', 'wb')
		pickle.dump(out_data, fil)
		fil.close()
		print 'Saved Cl datavector to disc.'

	def output_covariance_matrix(self, block):

		import pyfits

		cov = []
		order = []
		for spect1 in self.types:
			for spect2 in self.types:
				exec 'cov += [self.cov_%s%s]' %(spect1[0], spect2[0])
				order += ['%s%s'%(spect1[0], spect2[0])]

		filename = 'cl_covariance_matrix.fits'
		try: 
			outfile = pyfits.open(filename, mode='update')
		except:
			hdu=pyfits.PrimaryHDU( self.lbin_edges )
			hdulist = pyfits.HDUList([hdu])
			hdulist.writeto(filename)
			hdulist.close()
			outfile = pyfits.open(filename, mode='update')

		
			hdu = pyfits.PrimaryHDU(np.array(cov))
			hdu.header['omega_m'] = block['cosmological_parameters', 'omega_m']
			hdu.header['omega_de'] = block['cosmological_parameters', 'omega_de']
			hdu.header['h0'] = block['cosmological_parameters', 'h0'] 
			hdu.header['w0'] = block['cosmological_parameters', 'w0']
			hdu.header['wa'] = block['cosmological_parameters', 'wa']
			hdu.header['sigma_8'] = block['cosmological_parameters', 'sigma_8']

			outfile.append(hdu)
			oufile.flush()

#		cospar = {'omega_m' : block['cosmological_parameters', 'omega_m'], 
#			  'omega_de': block['cosmological_parameters', 'omega_de'],
#			  'h' : block['cosmological_parameters', 'h0'], 
#			  'w0' : block['cosmological_parameters', 'w0'],
#			  'wa' : block['cosmological_parameters', 'wa'],
#			  'sigma_8' : block['cosmological_parameters', 'sigma_8']}
 #
#		cov = { 'cov_eeee' : self.cov_eeee, 
#			'cov_eene' : self.cov_eene, 
#			'cov_eenn' : self.cov_eenn,
#			'cov_neee' : self.cov_neee,
#			'cov_nene' : self.cov_nene, 
#			'cov_nenn' : self.cov_nenn, 
#			'cov_nnee' : self.cov_nnee, 
#			'cov_nnne' : self.cov_nnne, 
#			'cov_nnnn' : self.cov_nnnn,
#			'l_bin_edges' : self.lbin_edges, 
#			'z_bin_edges' : self.zbin_edges}
		
#		import cPickle as pickle
#		fil = open(filename, 'wa')
#		pickle.dump(cov, fil)
#		fil.close()
		print 'Saved covariance matrix to disc.'
