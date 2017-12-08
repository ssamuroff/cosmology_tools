import numpy as np
import scipy as sp
import scipy.spatial as sps
import pyfits
from cosmosis.datablock import names, option_section
import os

import matplotlib.colors
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['font.size']=14
matplotlib.rcParams['legend.fontsize']=14
matplotlib.rcParams['xtick.major.size'] = 10.0
matplotlib.rcParams['ytick.major.size'] = 10.0

colours=['midnightblue','forestgreen','pink', 'purple', "lightblue"]
constants={"c" : 299792458.0}

def run_cosmosis(ini, values, extra_args="", sampler="test", outdir="cosmosis_output"):
	comm = "cosmosis %s -p pipeline.values=%s runtime.sampler=%s %s.save_dir=%s %s"%(ini,values,sampler,sampler,outdir,extra_args)
	print comm
	os.system(comm)

def postprocess(filename, outdir, contours=False, extra_flags="", factor_kde=1.6):
    print "Will postprocess chain : %s"%filename
    if not contours: contours = "--no-2d"
    else: contours = ""
    comm = "postprocess -o %s %s %s --factor-kde %2.2f %s"%(outdir, contours, extra_flags, factor_kde, filename)
    print "executing command : %s"%comm
    os.system(comm)

def get_double_or_array(options, param_name, section=option_section, default=0.0):
	if (section, param_name) not in options.keys():
		print "Found no parameter %s in section %s. Using default value ( %f ) "%(param_name, section, default)
		return default

	else:
		param = options[section, param_name]
		return param

def get_cosmological_parameters(block):
	try: omega_de = block['cosmological_parameters', 'omega_de']
	except:
		omega_k = block['cosmological_parameters', 'omega_k']
		omega_de = 1.0 - block['cosmological_parameters', 'omega_m'] - omega_k
	cospar={'omega_m': block['cosmological_parameters', 'omega_m'],
		'omega_de': omega_de,
		'omega_b': block['cosmological_parameters', 'omega_b'],
		'h': block['cosmological_parameters', 'h0'],
		'sigma_8': block['cosmological_parameters', 'sigma_8'],
		'n_s': block['cosmological_parameters', 'n_s'],
		'w0': block['cosmological_parameters', 'w'],
		'wa': block['cosmological_parameters', 'wa']
		}

	return cospar

def get_bias_parameters(block, sec="bias_parameters"):
	return {'b1': block[sec, 'bias_1'], 'b2': block[sec, 'bias_2'], 'b3': block[sec, 'bias_3']}

def get_galaxy_bias_parameters(block, sec="bias_parameters"):
	return {'b_1': block[sec, 'b_1'], 'b_2': block[sec, 'b_2'], 'b_3': block[sec, 'b_3']}

def get_ia_parameters(block):
	ia="intrinsic_alignment_parameters"
	return {'A_II': block[ia, 'A_II'], 'eta_II': block[ia, 'alpha_ii'], 'A_GI': block[ia, 'A_GI'], 'eta_GI': block[ia, 'alpha_GI']}

	
