import numpy as np
import scipy as sp
import pyfits
import os, pdb

class fisher:
	def __init__(self, fil):
		self.F = np.loadtxt(fil)
		self.param = []
		self.par = {}
		f = open(fil)
                varpar = f.readline().split("\t")
                varied = [i.replace("\n", "").replace('#', '') for i in varpar]
                print varied
		# Read in values config section
		val = f.read().split('VALUES_INI')[1]
		f.close()

		# Separate the lines
		val = val.replace('## ', '').split('\n')


		# Extract the varied parameters
		for c in val:
                        if '[' in c:
                            sec=c.replace('[', '').replace(']', '')
                            print "section: %s"%sec
			if ('=' in c) and (c.split(';')[0].count('.')>1):
				param_name = c.split('=')[0]. replace(' ', '')
				if '%s--%s'%(sec,param_name) not in varied:
                                    continue
				#if param_name in self.par.keys():
					#param_name+='_2'
				if "bias" in param_name:
					param_name=sec+"--"+param_name
                              

				self.param += [param_name]
				self.par[param_name] = {}
					
				param_range = c.split('=')[1].split(';')[0]
				i = 0
				for p in param_range.split(' '):
					try: 
						num = float(p)
						if i==0:
							self.par[param_name]['min'] = num
						elif i==1:
							self.par[param_name]['centre'] = num
						elif i==2:
							self.par[param_name]['max'] = num
						i=i+1
					except:
						continue
				
		print 'Found %d x %d matrix' %self.F.shape
		print 'Free parameters:', self.param 

	def remove_priors(self, param, sigma):
		"""Subtract off an existing prior from the Fisher matrix"""
		
		for i, p in enumerate(self.param):
			if p in param:
				j = np.argwhere(np.array(param)==p)[0,0]
				sub = 1.0 / sigma[j] / sigma[j]
				self.F[i,i] = self.F[i,i] - sub
				print "Removing prior on %s"%p

		print "Done"



	def reset_priors(self, priors = 'none', case = None):

		if os.path.exists(priors):
			print 'Reading priors from %s'%priors
			f = open(priors)
			lines = f.read()
			f.close()
			# Separate the lines
			lines = lines.replace('## ', '').split('\n')
		
			# Extract the varied parameters
			for l in lines:
				if ('=' in l) and (l.split(';')[0].count('.')>1):
					par = l.split('=')[0].replace(' ', '').lower()
					info = l.split('=')[1].split(' ')
					j = 0
					for i in info:
						#import pdb ; pdb.set_trace()
						if i=='':
							continue
						if j==0:
							prior_type= str(i)
							if prior_type.lower() != 'gaussian':
								print "Warning: specified prior type for %s is not 'gaussian'."%par
								print "The Fisher analysis can only use Gaussian priors."
							j+=1
						elif j==1:
							self.par[par]['prior_centre'] = float(i)
							j+=1
							print "Setting prior on %s, centre %f"%(par,float(i))
						elif j==2:
							self.par[par]['prior_sigma'] = float(i)

		elif case!=None:
			print 'Using inbuilt scenario : %s'%case

			# All parameters known exactly
			if case == 'delta':
				prior = 1.0e-35
				for p in self.param:
					self.par[p]['prior_sigma'] = prior

			# All parameters known exactly
			if case == 'pangloss':
				prior = 0.000001
				for p in self.param:
					self.par[p]['prior_sigma'] = prior

			#if case == 'optimistic':

			#if case == 'conservative':
		else:
			print 'No recognised option given. Please specify either a case (case = ) or a file name (priors = )'

	def get_priors(self, ignore=None):
		prior = []
		# Finally fit the priors together into a form that can be applied to 
		# the Fisher matrix
		ignore = np.atleast_1d(ignore)
		ig=False
		for p in self.param:
			if 'prior_sigma' in self.par[p].keys():
				if len(ignore)==0:
					ig=True
				elif p in ignore:
					prior+=[0.0]
				# This should read okay if a prior has been specified for this parameter.

				if (p not in ignore):
					prior += [self.par[p]['prior_sigma']**-2]
			else:
				# Otherwise set the corresponding matrix element to 0
				# Equivalent to a prior of infinite width (adds no new information)
				prior += [0.0]

		self.prior_information_matrix = np.diag(np.array(prior))
		assert self.prior_information_matrix.shape == self.F.shape

	def get_errorbar(self, parameter=None, case=None, priors='none', silent=False):
		if not parameter:
			raise ValueError('No valid parameter name selected.')

		if case!=None or priors!='none':
			self.reset_priors(case, priors)
			self.get_priors()
		else:
			if not silent:
				print 'Not resetting priors.'
				
		if not silent:
			print np.diag(self.prior_information_matrix)
		cov = np.linalg.inv(self.F + self.prior_information_matrix)

		# Get the covariance matrix element corresponding to sigma_8
		if isinstance(parameter, list):
			sel=[]
			for p in xrange(len(parameter)):
				sel.append("np.array(self.param)=='%s'"%parameter[p])
			i = np.argwhere(sel)
			x,y = i, i.T
			delta = np.linalg.det(cov[x,y])**0.5
		else:
			sel = "sel=np.array(self.param)=='%s'"%parameter
			exec sel
			i = np.argwhere(sel)
			delta = cov[i,i]**0.5
			delta = delta[0][0]

		if not silent: 
			print 'statistical errorbar = %f' %delta
		return delta

def photoz_prior_width_plot(fisher, dimension=1, lens=True, shear=True, show=False, sampling = None, bias_name='bias', bias_name2 = 's_z', parameter='sigma8_input'):
	nzbin = str(list(np.unique(fisher.param))).count(bias_name)

	print "%d z bins"%nzbin

	if "%s_1"%bias_name not in fisher.par.keys():
		print "Error: Fisher matrix does not have variable '%s'"%bias_name

	else:
		print 'Using bias type: %s'%bias_name
	if dimension==1:
		x = np.hstack((np.logspace(-60,-10,50), np.linspace(0.0001, 0.1, 100)))
                if bias_name == "s_z":
			x = np.linspace(0.00001, 0.2, 800)
		if sampling!=None:
			x = sampling
		delta_sigma_8 = np.zeros_like(x)
		for i, sig in enumerate(x):
			for bin in xrange(nzbin):
				#import pdb ; pdb.set_trace()
				# Set the photo-z prior width
				if shear:
				    fisher.par['%s_%d'%(bias_name, bin+1)]['prior_sigma'] = sig
				    if i==0 :
				    	print "Marginalising over %s"%bias_name
				if lens and ('%s_%d_2'%(bias_name, bin+1) in fisher.par.keys()):
					fisher.par['%s_%d_2'%(bias_name, bin+1)]['prior_sigma'] = sig
					if i==0:
						print 'Also marginalising over %s_%d_2'%(bias_name, bin+1)
				elif (not lens) and ('%s_%d_2'%(bias_name, bin+1) in fisher.par.keys()):
					fisher.par['%s_%d_2'%(bias_name, bin+1)]['prior_sigma'] = 1.0e-30
					if i==0:
						print 'Setting zero-width prior on lens bias %s_%d_2'%(bias_name, bin+1)
				

			fisher.get_priors(ignore=parameter)
			ds = fisher.get_errorbar(parameter = parameter, silent=np.invert(show))
			delta_sigma_8[i] = ds

	if dimension==2:
		x = np.linspace(0.0001, 0.06, 100)
		delta_sigma_8 = np.zeros((len(x), len(x)))
		for i, sig1 in enumerate(x):
			for j, sig2 in enumerate(x): 
				for bin in range(nzbin):
					# Set the photo-z prior width
					fisher.par['%s_%d'%(bias_name, bin+1)]['prior_sigma'] = sig1
					fisher.par['%s_%d'%(bias_name2, bin+1)]['prior_sigma'] = sig2  
				fisher.get_priors()
				if show: print "[%d,%d]"%(i,j), np.diag(fisher.prior_information_matrix)
				ds = fisher.get_errorbar(parameter = 'sigma8_input', silent=True)
				delta_sigma_8[i,j] = ds


	return np.array(x), np.array(delta_sigma_8)

def get_probe_comparison_plot(files_shear,files_ggl, files_pos, priors=None):

	Fs_shear=f.fisher(files_shear[0])
	Fsm_shear=f.fisher(files_shear[1])
	Ft_shear=f.fisher(files_shear[2])

	Fs_ggl=f.fisher(files_ggl[0])
	Fsm_ggl=f.fisher(files_ggl[1])
	Ft_ggl=f.fisher(files_ggl[2])

	Ft_pos=f.fisher(files_pos[2])
	Fs_pos=f.fisher(files_pos[0])
	Fsm_pos=f.fisher(files_pos[1])

	if priors!=None:
		Fs_shear.reset_priors(priors)
		Fs_ggl.reset_priors(priors)
		Fs_pos.reset_priors(priors)

		Fsm_pos.reset_priors(priors)
		Fsm_shear.reset_priors(priors)
		Fsm_ggl.reset_priors(priors)

		Ft_pos.reset_priors(priors)
		Ft_shear.reset_priors(priors)
		Ft_ggl.reset_priors(priors)

	ys=[]
	x,y = f.photoz_prior_width_plot(Fs_shear, dimension=1, bias_name='bias', show=False)
	ys+=[y]
	x,y = f.photoz_prior_width_plot(Fs_ggl, dimension=1, bias_name='bias', show=False)
	ys+=[y]
	x,y = f.photoz_prior_width_plot(Fs_pos, dimension=1, bias_name='bias', show=False)
	ys+=[y]

	ysm=[]
	x,y = f.photoz_prior_width_plot(Fsm_shear, dimension=1, bias_name='delta_z', show=False)
	ysm+=[y]
	x,y = f.photoz_prior_width_plot(Fsm_ggl, dimension=1, bias_name='s_z', show=False)
	ysm+=[y]
	x,y = f.photoz_prior_width_plot(Fsm_pos, dimension=1, bias_name='s_z', show=False)
	ysm+=[y]

	yt=[]
	x,y = f.photoz_prior_width_plot(Ft_shear, dimension=1, bias_name='t_z', show=False)
	yt+=[y]
	x,y = f.photoz_prior_width_plot(Ft_ggl, dimension=1, bias_name='t_z', show=False)
	yt+=[y]
	x,y = f.photoz_prior_width_plot(Ft_pos, dimension=1, bias_name='t_z', show=False)
	yt+=[y]

	import pylab as plt

	plt.subplot(311)
	plt.plot(x,ys[0]/0.82, 'm-', label='WL')
	plt.plot(x,ys[1]/0.82, 'g-', label='WL+ggl')
	plt.plot(x,ys[2]/0.82, 'b-', label='WL+ggl+LSS')
	plt.ylabel("fractional errorbar $\Delta \sigma_8/ \sigma_8$")
	plt.legend(loc='upper left')
	plt.xlim(0,0.1)
	plt.title("shift $\delta z$", loc='right')

	plt.subplot(312)
	plt.plot(x,ysm[0]/0.82, 'm-', label='WL')
	plt.plot(x,ysm[1]/0.82, 'g-', label='WL+ggl')
	plt.plot(x,ysm[2]/0.82, 'b-', label='WL+ggl+LSS')
	plt.ylabel("fractional errorbar $\Delta \sigma_8/ \sigma_8$")
	plt.title("smear $S_z$", loc='right')
	plt.xlim(0,0.1)

	plt.subplot(313)
	plt.plot(x,yt[0]/0.82, 'm-', label='WL')
	plt.plot(x,yt[1]/0.82, 'g-', label='WL+ggl')
	plt.plot(x,yt[2]/0.82, 'b-', label='WL+ggl+LSS')
	plt.title("tail $T_z$", loc='right')
	plt.xlim(0,0.1)
	plt.ylabel("fractional errorbar $\Delta \sigma_8/ \sigma_8$")
	plt.xlabel("prior width $\Delta p$")

	plt.show()




class fisherplots:
	def __init__(self, profile=None):
		if profile:
			self.profile(profile)

	def profile(self, profile):
		print 'Using plot profile %s'%profile
		import yaml
		config = yaml.load(open(profile))
		self.plots=[]
		self.sh= get_option( config['shear'], 'do_shear')
		if self.sh:
			self.plots+=['WL']
			print config['shear']['file']
		self.sh_ggl= get_option( config['shear_ggl'], 'do_shear_ggl')
		if self.sh_ggl:
			self.plots+=['WL+ggl']
			print config['shear_ggl']['file']
		self.sh_pos= get_option( config['shear_pos'], 'do_shear_pos')
		if self.sh_pos:
			self.plots+=['WL+LSS']
			print config['shear_pos']['file']
		self.sh_ggl_pos= get_option( config['shear_ggl_pos'], 'do_shear_ggl_pos')
		if self.sh_ggl_pos:
			self.plots+=['WL+ggl+LSS']
			print config['shear_ggl_pos']['file']
		
		print config['shear']['file']

		self.config = config

	def prior_plot(self, parameter='sigma8_input', lens=True, dimension=1, bias='bias', bias_2= 'none',photoz_only=False, plots=(1,1)):
		labels={'bias': '\delta z', 's_z': 'S_z', 't_z': 'T_z'}

		colours={'WL':'m-', 'WL+ggl':'g-', 'WL+LSS':'r-', 'WL+ggl+LSS':'b-'}

		import pylab as plt
		import sstools.fisher as f
		y=[]
		if self.sh:
			shear=f.fisher(self.config['shear']['file'])
			pr =  get_option(self.config['shear'], 'priors', default=None)
			if pr:
				shear.reset_priors(pr)

			if photoz_only:
				shear.reset_priors(case='delta')

                        print "shear-shear"
			x,ys = f.photoz_prior_width_plot(shear, lens=lens, bias_name=bias, parameter=parameter, dimension = dimension, bias_name2 = bias_2)
			y+=[ys]		
		if self.sh_ggl:
			sg=f.fisher(self.config['shear_ggl']['file'])
			pr =  get_option(self.config['shear_ggl'], 'priors', default=None)
			if pr:
				sg.reset_priors(pr)
			if photoz_only:
				sg.reset_priors(case='delta')
		
                        print "shear+ggl"
			x,ysg = f.photoz_prior_width_plot(sg, lens=lens, bias_name=bias, parameter=parameter, dimension = dimension, bias_name2 = bias_2)
			y+=[ysg]

		if self.sh_pos:
			sp=f.fisher(self.config['shear_pos']['file'])
			pr =  get_option(self.config['shear_pos'], 'priors', default=None)
			if pr:
				sp.reset_priors(pr)
			if photoz_only:
				sp.reset_priors(case='delta')

                        print "shear+pos"
			x,ysp = f.photoz_prior_width_plot(sp, lens=lens, bias_name=bias, parameter=parameter, dimension = dimension, bias_name2 = bias_2)
			y+=[ysp]
		if self.sh_ggl_pos:
			sgp=f.fisher(self.config['shear_ggl_pos']['file'])
			pr =  get_option(self.config['shear_ggl_pos'], 'priors', default=None)
			if pr:
				sgp.reset_priors(pr)
			if photoz_only:
				sgp.reset_priors(case='delta')

                        print "shear+ggl+pos"
			x,ysgp = f.photoz_prior_width_plot(sgp, lens=lens, bias_name=bias, parameter=parameter)
			y+=[ysgp]
		#import pdb ; pdb.set_trace()

		if parameter=='sigma8_input':
			f=0.82
		else:
			f=1.

		if dimension==2:
			return x, ys, ysg, ysp, ysgp


		# Statistical error relative to true value
                if plots[0]:
			for i, pl in enumerate(self.plots):
				plt.plot(x,y[i]/f, '%s'%colours[pl], label=pl )

			plt.legend(loc='upper left')
			plt.xlabel("prior width $\Delta %s$"%labels[bias])
			if parameter=='sigma8_input':
				plt.ylabel("fractional error $\Delta \sigma_8 / \sigma_8$")
			else:
				plt.ylabel("statistical error $\Delta %s$"%parameter)
			plt.show()
			plt.close()

		# degradation in errorbar relative to 0 prior width
		if plots[1]:
			for i, pl in enumerate(self.plots):
				plt.plot(x,y[i]/y[i][0], '%s'%colours[pl], label=pl)

			plt.legend(loc='upper left')
			plt.xlabel("prior width $\Delta %s$"%labels[bias])
			if parameter=='sigma8_input':
				plt.ylabel("degradation in statistical error $\Delta \sigma_8 / \Delta \sigma_{8,0}$")
			else:
				plt.ylabel("statistical error $\Delta %s / \Delta %s _0$"%(parameter,parameter))
			plt.show()

			plt.close()
		return 0

def get_option(config, option, default=''):
	if option in config.keys():
		return config[option]
	else:
		if default!='':
			return default
		else:
			print 'Cannot find option %s'%option
			return 0


