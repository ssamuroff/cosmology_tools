from cosmosis.datablock import names, option_section
import gaussian_cov as functions

import numpy as np
import matplotlib.pyplot as plt

def setup(options):
	shear = options[option_section, "shear"]
	position = options[option_section, "clustering"]
	survey = options[option_section, "survey"]
	
	config= {'shear': shear, 'clustering': position, 'survey': survey }
	return config

def execute(block, config):

	# Set up everything needed at the start
	Cl = functions.Cl_class(block, config)	
	Cl.load_Cls(block)

	# Then do the calculation
	Cl.calculate_all(block)

	#Cl.output_covariance_matrix(block)

	#import pdb ; pdb.set_trace()

	#Cl.plot_cov(mode=['e','e','e','e'])

	return 0

def cleanup(config):
	pass
