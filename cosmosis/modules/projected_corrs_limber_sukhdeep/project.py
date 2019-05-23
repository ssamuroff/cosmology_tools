from __future__ import print_function
from builtins import range
from cosmosis.datablock import names, option_section
import sys
import numpy as np
from scipy.interpolate import interp1d
from scipy.special import jn, jn_zeros
from scipy.interpolate import interp1d
from scipy.integrate import quad as scipy_int1d
from hankel_transform import *



def get_tapered_pk(spectrum_name, block, transform, kmax, kmin):
	P = block[spectrum_name, 'p_k_1_1']
	kh = block[spectrum_name, 'k_h']

	pk_taper = HT.taper(k=kh, pk=P, large_k_lower=5, large_k_upper=kmax,low_k_lower=kmin,low_k_upper=kmin*2)

	return pk_taper


def compute_c1_baseline():
    C1_M_sun = 5e-14  # h^-2 M_S^-1 Mpc^3
    M_sun = 1.9891e30  # kg
    Mpc_in_m = 3.0857e22  # meters
    C1_SI = C1_M_sun / M_sun * (Mpc_in_m)**3  # h^-2 kg^-1 m^3
    # rho_crit_0 = 3 H^2 / 8 pi G
    G = 6.67384e-11  # m^3 kg^-1 s^-2
    H = 100.
    H_SI = H * 1000.0 / Mpc_in_m  # h s^-1
    rho_crit_0 = 3 * H_SI**2 / (8 * np.pi * G)  # h^2 kg m^-3
    f = C1_SI * rho_crit_0
    return f

C1_rhoC = compute_c1_baseline()

def setup(options):
	rmax = options.get_double(option_section, 'rmax', default=200)
	rmin = options.get_double(option_section, 'rmin', default=0.1)
	kmax = options.get_double(option_section, 'kmax', default=75)
	kmin = options.get_double(option_section, 'kmin', default=1e-4)

	correlations = options[option_section, 'correlations'].split()

	do_bias = options.get_bool(option_section, 'add_bias', default=True)

	jn = []
	if 'wgp' in correlations:
		jn.append(2)
	if 'wgg' in correlations:
		jn.append(0)
	if 'wpp' in correlations:
		jn.append(0)
		jn.append(4)

	jn = np.unique(jn)

	print('Setting up Hankel transform...')
	print('this may take a little while')
	HT = hankel_transform(rmin=rmin, rmax=rmax, kmin=kmin, kmax=kmax, j_nu=[0,2,4], n_zeros=911000)

	return kmin, kmax, correlations, do_bias, HT


def execute(block, config):

	kmin, kmax, correlations, do_bias, HT = config

	# bias parameter, if needed 
	if do_bias:
		bg = block['bias_parameters', 'b0']
	else:
		bg = 1.


	# Hankel transforms, in turn
	# we've projected the power spectra along the z direction already,
	# so the P(k) here are 1D arrays 

	if 'wgg' in correlations:
		print('Processing wgg')

		P = get_tapered_pk('galaxy_power_projected', block, HT, kmax, kmin)
		r_gg,w_gg = HT.projected_correlation(k_pk=kh, pk=P*bg*bg, j_nu=0)

		block.put_double_array_1d('galaxy_w', 'w_rp_1_1', w_gg)
		block.put_double_array_1d('galaxy_w', 'r_p', r_gg)


	if 'wgp' in correlations:
		print('Processing wg+')
		P = get_tapered_pk('galaxy_intrinsic_power_projected', block, HT, kmax, kmin)
		r_gp,w_gp = HT.projected_correlation(k_pk=kh, pk=P*bg, j_nu=0)

		block.put_double_array_1d('galaxy_intrinsic_w', 'w_rp_1_1', w_gp)
		block.put_double_array_1d('galaxy_intrinsic_w', 'r_p', r_gp)
	if 'wpp' in correlations:
		P = get_tapered_pk('intrinsic_power_projected', block, HT, kmax, kmin)
		r_pp,w_pp = HT.projected_correlation(k_pk=kh, pk=P, j_nu=0)

		block.put_double_array_1d('intrinsic_w', 'w_rp_1_1', w_pp)
		block.put_double_array_1d('intrinsic_w', 'r_p', r_pp)


	return 0

