#Code by Donnacha Kirk
#Edited by Simon Samuroff 09/2015 

import numpy as np
from cosmosis.datablock import names as section_names
from cosmosis.datablock import option_section

def setup(options):
	dz = options.get_double(option_section, "dz", default=0.01)
	survey = options.get_string(option_section, "survey")
	return (dz, survey)

def gaussian(z,mu,sigma):
	return np.exp(-0.5*(z-mu)**2/sigma**2) / np.sqrt(2*np.pi) / sigma

def smail_distribution(z, alpha, beta, z0):
	return (z**alpha) * np.exp(-(z/z0)**beta)

def photometric_error(z, Nz, sigma_z, bias):
	nz = len(z)
	output = np.zeros((nz,nz))
	for i in xrange(nz):
		p = gaussian(z,z[i]-bias,sigma_z*(1+z[i]))
		#Could add f_cat or other terms here.
		#import pdb ; pdb.set_trace()
		output[:,i] = p * Nz[i]
	return output

def find_bins(z, nz_true, nbin):
	nz_true = nz_true/nz_true.sum()*nbin
	cum = np.cumsum(nz_true)
	bin_edges = [0.0]
	for i in xrange(1,nbin):
		edge = np.interp(1.0*i,cum,z)
		bin_edges.append(edge)
	bin_edges.append(z.max())	
	return np.array(bin_edges)

def compute_bin_nz(z_prob_matrix, z, edges, ngal):
	NI = []
	nbin = len(edges)-1
	dz = z[1]-z[0]
	for low,high in zip(edges[:-1], edges[1:]):
		w = np.where((z>low) & (z<high))[0]
		# Sum over all possible ztrue
		# Equivalent to marginalising p(zphot|ztrue) wrt ztrue
		ni = z_prob_matrix[w,:].sum(axis=0)

		# Normalise the n(z) in each redshift bin to 1 over the redshift range
		# of the survey
		ni*= 1.0/(ni.sum()*dz)							#ngal/(nbin*ni.sum()*dz)
		assert(len(ni)==len(z))
		NI.append(ni)
	return NI
	
def compute_nz(alpha, beta, z0, z, nbin, sigma_z, ngal, bias):
	#Set up Smail distribution of z vector as the distribution of true redshifts of the galaxies, n(ztrue)
	nz_true = smail_distribution(z, alpha, beta, z0)

	# Multiply that by a Gaussian to get the probability distribution of the measured photo-z for each true redshift p(zphot|ztrue)
	# This gives a 2D probability distribution
	z_prob_matrix = photometric_error(z, nz_true, sigma_z, bias)
	edges = find_bins(z,nz_true,nbin)
	bin_nz = compute_bin_nz(z_prob_matrix, z, edges, ngal)
	return edges,bin_nz

def execute(block, config):
	dz, survey = config
	params = section_names.number_density_params
	alpha = block[params, "alpha"]
	beta = block[params, "beta"]
	z0 = block[params, "z0"]

	sigma_z = block.get_double(survey, "sigz")
	bias = block.get(survey, "photoz_bias")
	ngal = block.get_double(survey, "ngal")
	zmax = block.get_double(survey, "zmax")
	# Using block.get_int() doesn't work as nzbin seems to be a double...
	nbin = int(block[survey, "nzbin"])
	
	#Compute the redshift vector
	z = np.arange(0,zmax+dz/2,dz)
	
	#Run the main code for getting n(z) in bins
	edges,bins = compute_nz(alpha, beta, z0, z, nbin, sigma_z, ngal, bias)

	#Save the results
	block[survey,"nz"] = len(z)
	block[survey,"z"] = z

	#Loop through the bins
	for i,bin in enumerate(bins):
		#The bin numbering starts at 1
		b=i+1
		name = "bin_%d" % b
		#Save the bin edges as parameters
		block[survey,"edge_%d"%b] = edges[i]
		#And save the bin n(z) as a column
		block[survey, name] =  bin
	#Also save the upper limit to the top bin
	block[survey,"edge_%d"%(nbin+1)] = edges[-1]

	return 0
		
def cleanup(config):
	#nothing to do here!  We just include this 
	# for completeness
	return 0
