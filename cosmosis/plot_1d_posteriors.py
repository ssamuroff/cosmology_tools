import numpy as np
from tools.likelihoods import *
import tools.emcee as mc
import sys
import tools.arrays as arr
import pylab as plt
plt.style.use("y1a1")
plt.switch_backend("pdf")

sys.path.append('/home/samuroff/cosmosis/')
from cosmosis.postprocessing import plots
from cosmosis.plotting import kde

files = sys.argv[2:]
print 'Found %d file(s) to process'%len(files)

def get_col(name, chain):
	if name in chain.samples.dtype.names:
		print name
		return chain.samples[name]
	else:
		return None

from matplotlib import rcParams
rcParams['figure.figsize'] = 10, 14

parameter_labels = {
	          'intrinsic_alignment_parameters--a_gi':r'$A_\mathrm{GI}$',
	           'intrinsic_alignment_parameters--a':r'$A_\mathrm{IA}$',
	          'intrinsic_alignment_parameters--a_ii':r'$A_\mathrm{II}$',
	          'intrinsic_alignment_parameters--c1':r'$A_1$',
	          'intrinsic_alignment_parameters--c2':r'$A_2$',
	          'intrinsic_alignment_parameters--alpha_ii':r'$\eta_\mathrm{GI}$',
	          'intrinsic_alignment_parameters--alpha':r'$\eta_\mathrm{IA}$',
	          'intrinsic_alignment_parameters--alpha_gi':r'$\eta_\mathrm{II}$',
	          'intrinsic_alignment_parameters--alpha_1':r'$\eta_1$',
	          'intrinsic_alignment_parameters--alpha_2':r'$\eta_2$',
	          'cosmological_parameters--s8':'$S_8$',
	          'intrinsic_alignment_parameters--bias_ta':r'$b_{g,s}$',
	          'intrinsic_alignment_parameters--bias_tt':r'$b_{g,s}$',
	          'intrinsic_alignment_parameters--a1':r'$A^{(1)}$',
	          'intrinsic_alignment_parameters--a2':r'$A^{(2)}$',
	          'intrinsic_alignment_parameters--a3':r'$A^{(3)}$',
	          'intrinsic_alignment_parameters--a4':r'$A^{(4)}$',
	          'cosmological_parameters--omega_m':r'$\Omega_\mathrm{m}$',
	          'cosmological_parameters--omega_b':r'$\Omega_\mathrm{b}$',
	          'cosmological_parameters--n_s':r'$n_\mathrm{s}$',
	          'cosmological_parameters--a_s':r'$A_\mathrm{s}$',
	          'cosmological_parameters--h0':r'$h$',
	          'cosmological_parameters--w':r'$w_0$',
	          'cosmological_parameters--sigma_8':r'$\sigma_8$',
	          'cosmological_parameters--omnuh2':r'$\Omega_\nu h^2$',
	          'wl_photoz_errors--bias_1':r'$\Delta z^{(1)}_{s}$',
	          'wl_photoz_errors--bias_2':r'$\Delta z^{(2)}_{s}$',
	          'wl_photoz_errors--bias_3':r'$\Delta z^{(3)}_{s}$',
	          'wl_photoz_errors--bias_4':r'$\Delta z^{(4)}_{s}$',
	          'shear_calibration_parameters--m1':r'$m^{(1)}$',
	          'shear_calibration_parameters--m2':r'$m^{(2)}$',
	          'shear_calibration_parameters--m3':r'$m^{(3)}$',
	          'shear_calibration_parameters--m4':r'$m^{(4)}$',
	          'lens_photoz_errors--bias_1':r'$\Delta z^{(1)}_{l}$',
	          'lens_photoz_errors--bias_2':r'$\Delta z^{(2)}_{l}$',
	          'lens_photoz_errors--bias_3':r'$\Delta z^{(3)}_{l}$',
	          'lens_photoz_errors--bias_4':r'$\Delta z^{(4)}_{l}$',
	          'lens_photoz_errors--bias_5':r'$\Delta z^{(5)}_{l}$',
	          'bias_parameters--b_1':r'$b^{(1)}_{g,l}$',
	          'bias_parameters--b_2':r'$b^{(2)}_{g,l}$',
	          'bias_parameters--b_3':r'$b^{(3)}_{g,l}$',
	          'bias_parameters--b_4':r'$b^{(4)}_{g,l}$',
	          'bias_parameters--b_5':r'$b^{(5)}_{g,l}$'
	          }

colours = ['purple', 'purple', 'purple', 'royalblue', 'royalblue', 'royalblue', 'purple', 'purple', 'purple']
ls = ['--', ':', '-']*3
alpha = [0., 0., 0.3]*3

panel_numbers = {} 

npar=None



for j,path in enumerate(files):
	chain = mc.chain(path)
	chain.add_column("s8", values="sigma_8*((omega_m/0.3)**0.5)")
	vals = arr.add_col(np.array(chain.samples),'weight',chain.weight)

	names = chain.samples.dtype.names
		
	npar = len(names)
	nrows = 1 + int(npar/5)
	
	ipanel = 0

	for i,name in enumerate(names):
		label = parameter_labels[name]
		if label not in panel_numbers.keys():
			ipanel+=1
			panel_numbers[label] = ipanel
		else:
			ipanel = panel_numbers[label]

		plt.subplot(nrows, 5, ipanel)

		x = get_col(name, chain)
		if x is None:
			continue
		n, x_axis, like = smooth_likelihood(x, kde)
		like/=like.max()

		plt.plot(x_axis, like, color=colours[j], ls=ls[j])
		if (alpha[j]>0.0):
			plt.fill_between(x_axis, like, color=colours[j], alpha=alpha[j])

		plt.xlabel(label, fontsize=14)
		plt.yticks(visible=False)
		plt.xticks(fontsize=8, rotation=45)

		if name=='cosmological_parameters--s8':
			plt.xlim(0.4,1.4)


plt.subplots_adjust(wspace=0, bottom=0.05, top=0.99, hspace=0.4, left=0.05, right=0.97)
plt.savefig(sys.argv[1])
