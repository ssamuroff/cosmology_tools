from cosmosis.datablock import names, option_section
import generate_observable_spectra as functions

import numpy as np
import matplotlib.pyplot as plt

def setup(options):
	shear= options[option_section, "shear"]
	intrinsic_alignments= options[option_section, "intrinsic_alignments"]
	clustering= options[option_section, "clustering"]
	magnification= options[option_section, "magnification"]
	Nl = options[option_section, "Nl"]
	lmin = options[option_section, "lmin"]
	lmax = options[option_section, "lmax"]
	noise = options[option_section, "noise"]
	
	opt= {'shear': shear, 'intrinsic_alignments': intrinsic_alignments, 'clustering': clustering, 'magnification': magnification, 'noise': noise, 'Nl': Nl, 'lmin': lmin, 'lmax': lmax }
	return opt

def execute(block, config):
	
	# Survey parameters
	Nl = config['Nl']
	lmax = config['lmax']
	lmin = config['lmin']

	# Defines whether to add noise 
	noise = config['noise']
	
	Cl = functions.Cl_class()
	Cl.load_and_generate_observable_cls(block, names, config)
	Cl.save_cls(block, config)

	return 0

def cleanup(config):
	pass
