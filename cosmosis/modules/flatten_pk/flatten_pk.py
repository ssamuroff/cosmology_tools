from __future__ import print_function
from builtins import range
from cosmosis.datablock import names, option_section
import sys
import numpy as np
import scipy.interpolate as spi
import scipy.integrate as sint
import scipy.special as sps

def setup(options):

    power_spectrum_name = options.get_string(option_section, "pk_name")
    sample_a = options.get_string(option_section, "sample_a", default="source")
    sample_b = options.get_string(option_section, "sample_b", default="lens")
    window_function = options.get_bool(option_section, "window_function", default=False)
    if not window_function:
    	redshift = options[option_section, "redshift"]
    	if isinstance(redshift, float):
    		redshift = [redshift]
    else:
    	redshift = None
    add_bias = options.get_bool(option_section, "add_bias", default=False)
    add_rsd = options.get_bool(option_section, "add_rsd", default=False)
    add_ia = options.get_bool(option_section, "add_intrinsic_alignments", default=False)

    return power_spectrum_name, redshift, add_bias, add_ia, add_rsd, window_function, sample_a, sample_b

def apply_ia(block, pkname, redshift, pk, add_ia):
	if (not add_ia) or (pkname=='galaxy_power'):
		return pk

	zname = '%3.3f'%(redshift)
	zname = zname.replace('.','_')
	A = block['intrinsic_alignment_parameters', 'a_%s'%zname]

	#print('Applying bg=%3.3f to %s'%(bg, pkname))

	if (pkname=='intrinsic_power'):
		return A * A * pk
	elif (pkname=='galaxy_intrinsic_power') or (pkname=='matter_intrinsic_power'):
		return A * pk
	else:
		raise ValueError('Unknown power spectrum type: %s'%pkname)

def apply_bias(block, pkname, redshift, pk, add_bias):
	if (not add_bias) or (pkname=='intrinsic_power'):
		return pk

	if not (redshift==None):
		zname = '%3.3f'%(redshift)
		zname = zname.replace('.','_')
		bg = block['bias_parameters', 'b_%s'%zname]

	else:
		bg = block['bias_parameters', 'b0']

	#print('Applying bg=%3.3f to %s'%(bg, pkname))

	if (pkname=='galaxy_power'):
		return bg * bg * pk
	elif (pkname=='galaxy_intrinsic_power'):
		return bg * pk
	else:
		raise ValueError('Unknown power spectrum type: %s'%pkname)

def apply_rsd(block, pkname, za, pk, add_rsd):
	if (not add_rsd) or (pkname=='intrinsic_power'):
		return pk

#	lnD = np.log(block['growth_parameters', 'd_z'])
#	z = block['growth_parameters', 'z']
#	lna = np.log(1./(1.+z))
#
#	lna0 = np.linspace(lna[0],lna[-1],2000)
#	lnD0 = spi.interp1d(lna, lnD, fill_value='extrapolate')(lna0)
#
#	da = lna0[0]-lna0[1]
#	fz = spi.interp1d(lna0, np.gradient(lnD0, da), fill_value='extrapolate')
#
	bg = block['bias_parameters', 'b0']
#
#	if (pkname=='galaxy_intrinsic_power'):
#		F = (1 + fz(za)/bg)
#	elif (pkname=='galaxy_power'):
#		F = (1 + fz(za)/bg) * (1 + fz(za)/bg) 


#	alpha = 0.55
#	fz = block['cosmological_parameters', 'omega_m'] * (1+za) * (1+za) * (1+za) 
#	fz = fz**alpha

	fz0 = block['growth_parameters', 'f_z']
	d_z = block['growth_parameters', 'd_z']
	z0 = block['growth_parameters', 'z']
	fz0 = (d_z / d_z[0]) * fz0
	interp = spi.interp1d(z0, fz0, fill_value='extrapolate')
	 
	fz = interp(za)
	betaD = interp(0.0)/bg

	if (pkname=='galaxy_intrinsic_power'):
		F = (1 + fz/bg)
	elif (pkname=='galaxy_power'):
		F = (1 + fz/bg) * (1 + fz/bg) 

	Fgrid,_ = np.meshgrid(F, pk[0])
	Fgrid = Fgrid.T




	k = block['matter_power_nl','k_h']

	Pk0_int = spi.interp1d(np.log(k), pk[0],fill_value='extrapolate')
	Pk0 = lambda k: Pk0_int(np.log(k))


	kPk = pk[0] * k
	Nk = len(kPk) ## number of bins

	kPk_fft = np.fft.fft(kPk) ## FFT
	klogstep = np.log(k[1]/k[0])

	kmin,kmax = k[[0,-1]]
	ifreq = np.fft.fftfreq(Nk)
	nu_arr = 1j*2*np.pi*np.fft.fftfreq(Nk, d=klogstep) ## the argument goes into k^nu_n, now nu_n is nu_arr
	cn = kPk_fft * kmin**(-nu_arr) ## the fourier coefficients

	kPk0 = np.array([np.sum(cn*k0**nu_arr)/Nk for k0 in k])

















	xi_gp = lambda kz,kp,rp,PI,z:np.cos(kz*PI)*kp**3/(kp**2+kz**2)*Pk0(kabs(kz,kp))*sps.jn(2,kp*rp)*(1.0+betaD*kz**2/(kp**2+kz**2))

	K = lambda kz,kp,z: np.cos(kz*PI)*kp**3/(kp**2+kz**2)*Pk0(kabs(kz,kp))*(1.0+betaD*kz**2/(kp**2+kz**2))

	xi_gp0 = lambda kz,kp,rp,PI,z:np.cos(kz*PI)*kp**3/(kp**2+kz**2)*Pk0(kabs(kz,kp))*sps.jn(2,kp*rp)

	nk = 501 ## similar answer between 500-5000 bins
	kabs = lambda kz, kp: np.sqrt(kz**2+kp**2)

	ik = np.logspace(-4, 2,nk)
	dk = ik[1:]-ik[:-1]
	ikc = 0.5*(ik[1:]+ik[:-1])
	kz, kp = np.array(np.meshgrid(ikc,ikc)).reshape(2,-1)
	dkz, dkp = np.array(np.meshgrid(dk,dk)).reshape(2,-1)

	def xi_int0 (rp, PI, z=0):
		xi_arr = xi_gp0(kz,kp,rp,PI,z)
		out = np.sum(xi_arr*dkz*dkp)
		return out

	def xi_int (rp, PI, z=0):
		xi_arr = xi_gp(kz,kp,rp,PI,z)
		out = np.sum(xi_arr*dkz*dkp)
		return out

	PI_arr = np.logspace(-0.5, 1, 50)
	rp_arr = np.logspace(-1,np.log10(200),25) #[1.0, 2.0, 5.0, 10.0, 20.0, 40.0]

	#import pdb ; pdb.set_trace()

	brute_arr0 = np.array([[xi_int0(irp, iPI) for iPI in PI_arr] for irp in rp_arr])
	brute_arr  = np.array([[xi_int(irp, iPI) for iPI in PI_arr] for irp in rp_arr])





	import pdb ; pdb.set_trace()

	#betaD = lambda z: interp(z)/bg
	#Fcoeff = lambda rp, PI, nu: 
	#rp**2/(rp**2+PI**2)**(2.0+nu/2.0)*cos(pi*nu/2.0)*gamma(3.0+nu)/(1-nu**2)


	#Fcoeff = lambda rp, PI, betaD: 
	#rp**2/(rp**2+PI**2)**(3.0+nu/2.0)*cos(pi*nu/2.0)*gamma(3.0+nu)*( (nu-3.0)*(rp**2+PI**2)+ (rp**2-(3+nu)*PI**2*betaD))/(nu-3.0)/(nu**2-1)

	return Fgrid * pk


def execute(block, config):
	power_spectrum_name, redshift, add_bias, add_ia, add_rsd, window_function, sample_a, sample_b = config

	z,k,pk = block.get_grid(power_spectrum_name, 'z', 'k_h', 'p_k')

	if len(pk[pk>0])==len(pk):
		logint = True
		lnpk = np.log(pk)
	else:
		print('Negative Pk values - will interpolate in linear space.')
		logint = False
		lnpk = pk

	interp = spi.interp2d(np.log(k), z, lnpk, kind='linear')

	# Simplest case: take the 2D grid P(k,z) and interpolate to the desired redshift
	# Effectively assuming the redshift distribution is a delta fn.
	# Store the linearised power spectrum to the block with the suffix identifying z
	# Only do this if window_function=F
	redshift = np.atleast_1d(redshift)
	if (not window_function):
		for z0 in redshift:
			pk_interpolated = interp(np.log(k), [z0])
			if logint:
				pk_interpolated = np.exp(pk_interpolated)

			pk_interpolated = apply_ia(block, power_spectrum_name, z0, pk_interpolated, add_ia)
			pk_interpolated = apply_bias(block, power_spectrum_name, z0, pk_interpolated, add_bias)

			block.put_double_array_1d(power_spectrum_name+'_%2.3f'%z0, 'p_k', pk_interpolated)
			block.put_double_array_1d(power_spectrum_name+'_%2.3f'%z0, 'k_h', k)
			#import pdb ; pdb.set_trace()

	# Slightly more complicated case: project P(k,z) along the redshift axis,
	# with a window function defined by a specific pair of n(z)
	# This is a bit more fiddly as it involves looping over bin pairs.
	else:

		# we'll need chi(z) and dchi/dz
		# (at the same z sampling as the n(z) )
		chi = block['distances','d_m']
		z = block['distances','z']
		interp_chi = spi.interp1d(z,chi)
		dz = z[1]-z[0]
		Dchi = np.gradient(chi,dz)
		interp_dchi = spi.interp1d(z,Dchi)

		# Now do the kernel calculation for each bin pair in turn
		# see, for example, eq 8 from Singh et al arXiv:1411.1755
		na = block['nz_%s'%sample_a, 'nbin']
		nb = block['nz_%s'%sample_b, 'nbin']
		zmin = 0.01

		# Do this only once per power spectrum
		output_section_name = str(power_spectrum_name+'_projected')
		block.put_double_array_1d(output_section_name, 'k_h', k)
		block.put_int(output_section_name, 'nbin_a', na)
		block.put_int(output_section_name, 'nbin_b', nb)

		#Now loop over bin pairs
		for i in range(na):
			nz_a = block['nz_%s'%sample_a, 'bin_%d'%(i+1)]
			za = block['nz_%s'%sample_a, 'z']
			for j in range(nb):
				nz_b = block['nz_%s'%sample_b, 'bin_%d'%(j+1)]
				zb = block['nz_%s'%sample_b, 'z']

				if len(za)!=len(zb):
					interp_nz = spi.interp1d(zb, nz_b)
					nz_b = interp_nz(za)
					#raise ValueError('Redshift sampling does not match!')

				x = interp_chi(za)
				dxdz = interp_dchi(za)

				X = nz_a * nz_b/x/x/dxdz
				interp_X = spi.interp1d(za, X)

				# Inner integral over redshift
				V,Verr = sint.quad(interp_X, zmin, za.max())
				W = nz_a*nz_b/x/x/dxdz/V

				# Now do the power spectrum integration
				W2d,_ = np.meshgrid(W,k)
				W2d[np.invert(np.isfinite(W2d))] = 1e-30

				pk_interpolated = interp(np.log(k), za)
				if logint:
					pk_interpolated = np.exp(pk_interpolated)

				pk_interpolated = apply_bias(block, power_spectrum_name, None, pk_interpolated, add_bias)
				pk_interpolated = apply_rsd(block, power_spectrum_name, za, pk_interpolated, add_rsd)

				# Outer integral with kernel calculated abov
				integrand = W2d.T * pk_interpolated
				Pw = sint.trapz(integrand,za,axis=0)

				block.put_double_array_1d(output_section_name, 'p_k_%d_%d'%(i+1,j+1), Pw)

	return 0





