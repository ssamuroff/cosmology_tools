import numpy as np
#from tools.iacorrs.power_spectra import *
import scipy
from scipy.interpolate import interp1d
from scipy.special import jn, jn_zeros
from scipy.interpolate import interp1d
from scipy.integrate import quad as scipy_int1d
from tools.iacorrs.hankel_transform import *
from astropy.cosmology import Planck15 #use Planck15 if you can
import astropy.units as u
import sys
import yaml

from matplotlib import rc
rc('text', usetex=False)

path = '/home/ssamurof/mbii/chains/example_output/'
#path = '/Users/hattifattener/coma/mbii/chains/example_output/'

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
#0.0134

def DZ_int(z=[0],cosmo=None,rtol=1.e-4,tol=1.e-5):
    #linear growth factor.. full integral.. eq 63 in Lahav and Suto
    Ez_func=cosmo.efunc
    def intf(z):
    	return (1.+z)/(cosmo.H(z=z).value**3)
    dz=np.zeros_like(z,dtype='float32')
    inf=np.inf
    j=0
    for i in z:
    	dz[j]+=cosmo.H(i).value*scipy_int1d(intf,i,inf,epsrel=rtol,epsabs=tol)[0]
    	j=j+1
    dz=dz*2.5*cosmo.Om0*cosmo.H0**2
    return dz/dz[0] #check for normalization

import argparse
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--config', type=str, action='store')
parser.add_argument('--redshift', type=float, action='store', default=-1.0)
parser.add_argument('--process', type=str, action='store', nargs='+', required=True)
parser.add_argument('--ind', type=int, action='store')
args = parser.parse_args()


settings = yaml.load(open(args.config,'rb'))

correlations = args.process
print('%d correlations to process: '%len(correlations), correlations)
nc = len(correlations)

# Decant into dictionaries
cosmo_fid = {}
for name in settings['cosmology'].keys():
	cosmo_fid[name] = settings['cosmology'][name]
	print(name, cosmo_fid[name])

pk_params=settings['general']
cosmo=Planck15.clone()
#cosmo_fid=dict({'h':cosmo.h,'Omb':cosmo.Ob0,'Omd':cosmo.Om0-cosmo.Ob0,'s8':0.817,'Om':cosmo.Om0,'As':2.12e-09,'mnu':cosmo.m_nu[-1].value,'Omk':cosmo.Ok0,'tau':0.06,'ns':0.965,'w':-1,'wa':0})

# Generate the matter power spectrum
# using CCL rather than class, just because this works on my laptop
if (args.redshift==-1.0):
	zlim = settings['general']['redshift']
else:
	zlim = args.redshift
print('Generating theory for z=%3.3f'%zlim)

#PS=Power_Spectra(cosmo_params=cosmo_fid,pk_params=pk_params)
z = [zlim,zlim+0.001]

#pk,kh =PS.class_pk(z=z)
#pk,kh = PS.ccl_pk(z=z)


pk = [np.loadtxt('%s/galaxy_power_%1.3f/p_k.txt'%(path,z[0]))]
kh = np.loadtxt('%s/galaxy_power_%1.3f/k_h.txt'%(path,z[0]))

#Setting up the Hankel Transform
#This part is slower. But only needs to be run once. 
#If you only need wgg, set j_nu=[0]. For wg+ (or \Delta\Sigma) use j_nu=[2]
jnu = []
if 'wgg' in correlations:
	jnu.append(0)
if 'wgp' in correlations:
	jnu.append(2)
if 'wpp' in correlations:
	jnu.append(4)

HT = hankel_transform(rmin=settings['general']['rmin'],rmax=settings['general']['rmax'],kmax=settings['general']['kmax'],j_nu=jnu,n_zeros=111000,kmin=settings['general']['kmin'])
#HT=hankel_transform(rmin=1,rmax=rmax,kmax=1,j_nu=[0,2],n_zeros=2800,kmin=1.e-2)#quick test... inaccurate

pk_taper=HT.taper(k=kh, pk=pk[0],large_k_lower=5, large_k_upper=settings['general']['kmax'],low_k_lower=settings['general']['kmin'],low_k_upper=settings['general']['kmin']*2)
#need atleast k=10 and k=1.e-3 to get decent wgg for 1-100 Mpc/h range.

# Get the coefficients
omega_m = settings['cosmology']['Om']
Dz = DZ_int(z=np.append([0],z),cosmo=cosmo)

#Dzc = np.loadtxt('/Users/hattifattener/coma/mbii/chains/example_output/growth_parameters/d_z.txt')
#zc = np.loadtxt('/Users/hattifattener/coma/mbii/chains/example_output/growth_parameters/z.txt')
#interpolator = interp1d(zc,Dzc)
#Dz = interpolator(np.append([0],z))
if 'overwrite' in settings['nuisance'].keys():
	path = settings['nuisance']['overwrite']
	print('Reading parameter values from %s'%path)
	settings['nuisance']['bias'], settings['nuisance']['a_ia'] = np.loadtxt(path)[args.ind], np.loadtxt(path)[args.ind+1]

wgg_f = settings['nuisance']['bias']**2
wgp_f = settings['nuisance']['bias']*settings['nuisance']['a_ia']*C1_rhoC*omega_m/Dz
wpp_f = (settings['nuisance']['a_ia']*C1_rhoC*omega_m/Dz)**2

def export(filename, x, y):
	out = np.vstack((x,y))
	np.savetxt(filename, out.T)
	print('Saved to %s'%filename)

# Now do the actual correlations
if 'wgg' in correlations:
	print('Processing wgg')
	#r_gg,wgg=HT.projected_correlation(k_pk=kh,pk=pk[0]*wgg_f,j_nu=0)
	r_gg,wgg_taper = HT.projected_correlation(k_pk=kh,pk=pk_taper*wgg_f,j_nu=0)
	interp_gg = interp1d(np.log10(r_gg), np.log10(wgg_taper), fill_value='extrapolate')
	R = np.array([ 0.11560133,  0.15448579,  0.20644969,  0.27589253,  0.36869363, 0.49270995,  0.65844125,  0.87991906,  1.17589468,  1.57142668, 2.10000253,  2.80637379,  3.75034493,  5.01183667,  6.69765243, 8.95052074, 11.9611793 , 15.98452361, 21.36118761, 28.5463832 ])
	export('wgg_%3.3f.txt'%zlim, R, 10**(interp_gg(np.log10(R))))


if 'wgp' in correlations:
	print('Processing wg+')
	r_gp,wgp = HT.projected_correlation(k_pk=kh,pk=pk_taper*wgp_f[1],j_nu=2)
	interp_gp = interp1d(np.log10(r_gp), np.log10(wgp), fill_value='extrapolate')
	R = np.array([ 0.13363669,  0.2386586 ,  0.42621476,  0.76116687,  1.35934994, 2.42763096,  4.33544879,  7.74257561, 13.82728292, 24.693818  ])
	export('wgp_%3.3f.txt'%zlim, R, 10**(interp_gp(np.log10(R))))

if 'wpp' in correlations:
	print('Processing w++')
	r_pp,wpp4 =  HT.projected_correlation(k_pk=kh,pk=pk_taper*wpp_f[1],j_nu=4)
	r_pp0,wpp0 = HT.projected_correlation(k_pk=kh,pk=pk_taper*wpp_f[1],j_nu=0)
	wpp0_intp = interp1d(r_pp0,wpp0,bounds_error=False,fill_value='extrapolate')
	wpp = wpp4 + wpp0_intp(r_pp) 
	interp_pp = interp1d(np.log10(r_pp), np.log10(wpp), fill_value='extrapolate')
	R = np.array([ 0.13363669,  0.2386586 ,  0.42621476,  0.76116687,  1.35934994, 2.42763096,  4.33544879,  7.74257561, 13.82728292, 24.693818  ])
	export('wpp_%3.3f.txt'%zlim, R, 10**(interp_pp(np.log10(R))))


print('Done all.')
