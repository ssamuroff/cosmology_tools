import numpy as np
import scipy as sp
import pyfits
import os, pdb
from sstools import fisher as fi
from samplers import sampler

import matplotlib.colors
import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['font.size']=14
matplotlib.rcParams['legend.fontsize']=14
matplotlib.rcParams['xtick.major.size'] = 10.0
matplotlib.rcParams['ytick.major.size'] = 10.0

colours=['midnightblue','forestgreen','pink', 'purple', "lightblue"]

param_names={"sigma8_input" : ["cosmological_parameters--sigma8_input", "\sigma_8"], "omega_m": ["cosmological_parameters--omega_m", "\Omega_m"],  "w" : ["cosmological_parameters--w", "w_0"], "wa" : ["cosmological_parameters--wa", "w_a"]}

fiducial={"sigma8_input" : 0.82, "omega_m":0.312, "w" : -1.0, "wa" : 0.0}
class grid(sampler):
	def __init__(self, fil):
		sampler.__init__(self,fil)
                self.extract_and_reshape()
                self.file = fil


        def extract_and_reshape(self):
            """Extract the posterior sample grid and reshape into a 2d matrix"""
            self.posterior = self.F.T[-1]
            self.nsamp = len(np.unique( self.F.T[0] ))
            self.npar = len(self.param)
            shape = tuple(self.npar*[self.nsamp])
	    self.posterior = self.posterior.reshape(shape)
		
	    # Do the same for the sample coordinates in each dimension
	    for i, p in enumerate(self.param):
                self.par[p]['samples'] = self.F.T[i].reshape(shape)
                self.par[p]['marginalised'] = False


	def marginalise(self, parameter, *args, **kwargs):
		"""Marginalise over the specified dimension in parameter space."""

		mode = kwargs.get('mode', 'internal')
		# index to marginalise over
		m = np.argwhere(np.array(self.param)==parameter)[0,0]

		import scipy.misc as spm
		samp = np.unique(self.par[parameter]['samples'])
		dp = samp[1]-samp[0]
		# Sum PDF along one axis and output the result in the required manner
		if mode=='internal':
			self.posterior = spm.logsumexp(self.posterior, axis=m)
			print self.posterior.shape

			# Clean up the parameter array
			self.param= np.delete(np.array(self.param), m)
			self.par[parameter]['marginalised'] = True
			print 'Marginalised over %s.'%parameter

			return 0

		elif mode=='return':
			posterior = args[0]
			posterior = spm.logsumexp(posterior, axis=m)
			param = args[1]
			param = np.delete(np.array(param), m)
			print 'Marginalised over %s.'%parameter

			return posterior, param		

	def get_1sigma_errorbar(self, parameter, frac=0.68):
		"""Repeatedly integrate a marginalised posterior distribution to
		   evaluate the 1 sigma statistical errorbar. """

		# Extract the sample points
		x = np.unique(self.par[parameter]['samples'])
		p_of_x = self.posterior
		# Resample from the specified parameter range using a higher resolution
		xf = np.linspace(x.min(), x.max(), len(x)*20)
		dx = xf[1]-xf[0]
		pf = np.interp(xf, x, p_of_x)
		pf = np.exp(pf)

		prob = sum(pf*dx)

		# Find the PDF peak
		i0 = np.argwhere(pf == pf.max())[0,0]
		x0 = xf[i0]
		print 'Integrating around %s = %f' %(parameter, x0)

		f = 0.0
		sigma = 0.0
		# Repeatedly integrate the PDF, extending the integration limits slightly
		# on each iteration
		for i in xrange(1, len(xf)):
			imin = i0-i
			imax = i0+i
			integrand = pf[imin:imax]
			f = sum(integrand * dx) / prob
			# If these bounds contain 68% of the total probability, save the
			# parameter range and break the loop
			if f>=frac:
				sigma = (xf[i0+i] - xf[i0-i])/2.
				print 'Estimated 1 sigma uncertainty range = %f after %d iterations' %(sigma,i)
				break

		return x0, sigma

class multigrid:
    def __init__(self, grids):
        for i, g in enumerate(grids):
            print "loading %s"%g
            setattr(self, "grid%d"%i, grid(g))
        self.names=grids

    def get_means(self, par):
        self.means=[]
        for i, g in enumerate(self.names):
            self.means+=[np.trapz(getattr(self,"grid%d"%i).par[par][ "samples"]*getattr(self,"grid%d"%i).posterior)/np.trapz(getattr(self,"grid%d"%i).posterior)]

    def savetotxt(self, filename, biases, par):
        self.get_means(par)
        out = np.array([biases, self.means]).T
        np.savetxt(filename,out) 

def pp_to_bias(loc, fil, method="best_fit"):
    lines =open("%s/means.txt"%loc).read()
    lines = lines.split("\n")

    b=[]
    p=[]
    p_std=[]
    for l in lines[1:]:
        if not l:
            continue
        if "#" in l:
            print l, l.split("zbias")[1].replace(".txt","")
            b+=[float(l.split("zbias")[1].replace(".txt",""))]
        else:
            
            r = l.split(" ")
            print r
            p+=[r[3]]
            p_std+=[r[6]]

    if method=="best_fit":
    	b=[]
    	p=[]
        lines =open("%s/best_fit.txt"%loc).read()
        lines = lines.split("\n")
        for l in lines[1:]:
            if not l:
                continue
            if "#" in l:
                b+=[float(l.split("zbias")[1].replace(".txt",""))]
            elif l.split(" ")[0]=="cosmological_parameters--sigma8_input":
                r= l.split(" ")
                p+=[r[-1]]

    order = np.argsort(np.array(b).astype(float))

    out=np.array([np.array(b).astype(float)[order], np.array(p).astype(float)[order], np.array(p_std).astype(float)[order]])
    print out
    np.savetxt(fil, out.T)

def pp_to_error(loc, fil, save=True):

	lines =open("%s/means.txt"%loc).read()
	lines = lines.split("\n")

	b=[]
	error={}
	for i, line in enumerate(lines[1:]):
		if ".txt" in line: 
			b+=[float(line.split("zprior")[1].split(".txt")[0])]
		elif "--" in line:
			param=line.split()[0]
			if param not in error.keys():
				error[param]=[float(line.split()[2])]
			else: error[param]+=[float(line.split()[2])]

	order = np.argsort(np.array(b).astype(float))
	out = np.array([np.array(b).astype(float)[order]])

	for n in error.keys():
		out=np.vstack((out, np.array(error[n]).astype(float)))
		error[n]=np.array(error[n])

	if save:
		np.savetxt(fil, out.T)
		f=open(fil, "a")
		f.write(lines[0].replace("#", "#zprior ")+"\n")
		f.close()

	return b, error


def pp_to_error_2d(loc, fil, save=True):

	lines =open("%s/ellipse_areas.txt"%loc).read()
	lines = lines.split("\n")

	b=[]
	error={}
	for i, line in enumerate(lines[1:]):
		if ".txt" in line: 
			b+=[float(line.split("zprior")[1].split(".txt")[0])]
		elif "--" in line:
			param1=line.split()[0]
			param2=line.split()[1]
			if (param1, param2) not in error.keys():
				error[(param1, param2)]=[float(line.split()[2])]
			else:
				error[(param1, param2)]+=[float(line.split()[2])]

	order = np.argsort(np.array(b).astype(float))
	out = np.array([np.array(b).astype(float)[order]])

	for n in error.keys():
		out=np.vstack((out, np.array(error[n]).astype(float)))
		error[n]=np.array(error[n])

	if save:
		np.savetxt(fil, out.T)
		f=open(fil, "a")
		f.write(lines[0].replace("#", "#zprior ")+"\n")
		f.close()

	return b, error

class fig5:
	def __init__(self, postprocess_dir, sec=None, shear=True, sp=True, sg=False, sgp=True, gp=True, from_fisher=False, prior=None, param="sigma8_input", do2d=False):
		self.dir=postprocess_dir
		if not from_fisher:
		    if shear:
		    	self.shear_shear = pp_to_error("%s/shear"%postprocess_dir, None, save=False)
		    if sp:
		    	self.sp = pp_to_error("%s/shear+pos"%postprocess_dir, None, save=False)		
		    if sg:
		    	self.sg = pp_to_error("%s/shear+ggl"%postprocess_dir, None, save=False)
		    if sgp:
		    	self.sgp = pp_to_error("%s/shear+ggl+pos"%postprocess_dir, None, save=False)
		    if gp:
		    	self.gp = pp_to_error("%s/ggl+pos"%postprocess_dir, None, save=False)

		if from_fisher:
			self.process_fisher(from_fisher, prior, shear, sp, sgp, gp, sg, sec, do2d, param )

	def process_fisher(self, filenames, prior, shear=True, sp=True, sgp=True, gp=True, sg=False, sec=None, do2d=False, param="sigma8_input"):
		for filename in filenames:
			f = fi.fisher(filename)
			if not sec:
				sec = "nz_shear_errors--bias"

			nparam = str(list(np.unique(f.param))).count(sec)
			sections=[]
			for i in xrange(nparam):
				sections.append(sec+"_%d"%(i+1))
			f.remove_priors(sections, np.array(nparam*[prior]))

			if param not in f.param:
				param="omega_m"
			x, y = fi.photoz_prior_width_plot(f, bias_name=sec, parameter=param)

			err = {param_names[param][0]:y}
			

			if "shear+ggl+pos" in filename:
				self.sgp = [x, err]
			if "shear+pos" in filename:
				self.sp = [x, err]
			if "ggl+pos" in filename:
				self.gp = [x, err]
			if ("shear+ggl" in filename) and ("pos" not in filename):
				self.sg = [x, err]
			if ("ggl" not in filename) and ("pos" not in filename):
				self.shear_shear = [x, err]

			if do2d:
				err1, err2 = {}, {}
				for p1 in f.param:
					for p2 in f.param:
						if p1==p2:
							continue
						else:
							x, err1[(p1,p2)] = fi.photoz_prior_width_plot(f, bias_name=sec, parameter=[p1,p2])
							x, err2[(p2,p1)] = fi.photoz_prior_width_plot(f, bias_name=sec, parameter=[p1,p2])
						if "shear+ggl+pos" in filename:
							self.sgp_2d = [x, err1]
							self.sgp_2d = [x, err2]
						if "shear+pos" in filename:
							self.sp_2d = [x, err1]
							self.sp_2d = [x, err2]
						if "ggl+pos" in filename:
							self.gp_2d = [x, err1]
							self.gp_2d = [x, err2]
						if ("shear+ggl" in filename) and ("pos" not in filename):
							self.sg_2d = [x, err1]
							self.sg_2d = [x, err2]
						if ("ggl" not in filename) and ("pos" not in filename):
							self.shear_shear_2d = [x, err1]
							self.shear_shear_2d = [x, err2]


	def contour_areas(self, shear=True, sp=True, sgp=True, gp=True):
		if shear:
			self.shear_shear_2d = pp_to_error_2d("%s/shear"%self.dir, None, save=False)
		if sp:
			self.sp_2d = pp_to_error_2d("%s/shear+pos"%self.dir, None, save=False)
		if sgp:
			self.sgp_2d = pp_to_error_2d("%s/shear+ggl+pos"%self.dir, None, save=False)
		if gp:
			self.gp_2d = pp_to_error_2d("%s/ggl+pos"%self.dir, None, save=False)

		print "Loaded countour areas from %s"%self.dir



	def make(self, loc=None, normalise=False, xlim_upper=0.06, ylim_lower=None, ylim_upper=None, s=True, sp=True, gp=True, sg=False, sgp=True, param="sigma8_input"):
		import pylab as plt

		if isinstance(param, list):
			name="nd_contour_area"
			for p in param:
				name+="--%s"%p
		else:
			name = param_names[param][0]

		if normalise:
			if s:
				plt.plot(self.shear_shear[0], self.shear_shear[1][name]/self.shear_shear[1][name][0], colours[0], label="WL", linewidth=2.0)
			if sp:
				plt.plot(self.sp[0], self.sp[1][name]/self.sp[1][name][0], colours[1], label="WL+LSS", linewidth=2.0)
			if gp:
				plt.plot(self.gp[0], self.gp[1][name]/self.gp[1][name][0], colours[2], label="GGL+LSS", linewidth=2.0)
			if sg:
				plt.plot(self.sg[0], self.sg[1][name]/self.sg[1][name][0], colours[3], label="WL+GGL", linewidth=2.0)
			if sgp:
				plt.plot(self.sgp[0], self.sgp[1][name]/self.sgp[1][name][0], colours[4], label="WL+GGL+LSS", linewidth=2.0)
			plt.ylabel("error degradation $\Delta %s/ \Delta %s (\Delta \delta z=0)$"%(param_names[param][1], param_names[param][1]))

		else:
			norm = abs(fiducial[param])
			if norm==0.0:
				norm=1.0
				ab=True
			else:
				ab=False
			if s:
				plt.plot(self.shear_shear[0], self.shear_shear[1][name]/norm, colours[0], label="WL", linewidth=2.0)
			if sp:
				plt.plot(self.sp[0], self.sp[1][name]/norm, colours[1], label="WL+LSS", linewidth=2.0)
			if gp:
				plt.plot(self.gp[0], self.gp[1][name]/norm, colours[2], label="GGL+LSS", linewidth=2.0)
			if sg:
				plt.plot(self.sg[0], self.sg[1][name]/norm, colours[3], label="WL+GGL", linewidth=2.0)
			if sgp:
				plt.plot(self.sgp[0], self.sgp[1][name]/norm, colours[4], label="WL+GGL+LSS", linewidth=2.0)
			if not ab: 
				plt.ylabel("fractional error $\Delta %s / %s$"%(param_names[param][1], param_names[param][1]))
			else:
				plt.ylabel("absolute error $\Delta %s$"%(param_names[param][1]))


		plt.xlabel("prior width $\Delta \delta z$")
		plt.legend(loc="upper left")
		plt.xlim(xmax=xlim_upper)
		plt.ylim(ymin=ylim_lower, ymax=ylim_upper)

		if not loc:
			plt.show()
		else:
			plt.savefig(loc)

		plt.close()



	def make2d(self, loc=None, normalise=False, xlim_upper=0.06, ylim_lower=None, ylim_upper=None, param=["sigma8_input", "omega_m"]):
		import pylab as plt

		p1 = param[0]
		p2 = param[1]

		if normalise:
			plt.plot(self.shear_shear_2d[0], self.shear_shear_2d[1][(p1,p2)]/self.shear_shear_2d[1][(p1,p2)][0], colours[0], label="WL", linewidth=2.0)
			plt.plot(self.sp_2d[0], self.sp_2d[1][(p1,p2)]/self.sp_2d[1][(p1,p2)][0], colours[1], label="WL+LSS", linewidth=2.0)
			plt.plot(self.gp_2d[0], self.gp_2d[1][(p1,p2)]/self.gp_2d[1][(p1,p2)][0], colours[2], label="GGL+LSS", linewidth=2.0)
			plt.plot(self.sgp_2d[0], self.sgp_2d[1][(p1,p2)]/self.sgp_2d[1][(p1,p2)][0], colours[3], label="WL+GGL+LSS", linewidth=2.0)
			#pdb.set_trace()
			plt.ylabel("error degradation $FOM^{-1}_{%s %s} / FOM^{-1}_{%s %s} (\Delta \delta z=0)$"%(param_names[param[0]][1], param_names[param[0]][1], param_names[param[1]][1], param_names[param[1]][1]))
		else:
			plt.plot(self.shear_shear_2d[0], self.shear_shear_2d[1][(p1,p2)], colours[0], label="WL", linewidth=2.0)
			plt.plot(self.sp_2d[0], self.sp_2d[1][(p1,p2)], colours[1], label="WL+LSS", linewidth=2.0)
			plt.plot(self.gp_2d[0], self.gp_2d[1][(p1,p2)], colours[2], label="GGL+LSS", linewidth=2.0)
			plt.plot(self.sgp_2d[0], self.sgp_2d[1][(p1,p2)], colours[3], label="WL+GGL+LSS", linewidth=2.0)
			#plt.yscale("log")
			plt.ylabel("$FOM^{-1}_{%s %s}$"%(param_names[p1][1], param_names[p2][1]))

		plt.xlabel("prior width $\Delta \delta z$")
		plt.legend(loc="upper left")
		plt.xlim(xmax=xlim_upper)
		plt.ylim(ymin=ylim_lower, ymax=ylim_upper)

		if not loc:
			plt.show()
		else:
			plt.savefig(loc)

		plt.close()








        

#def bias_plots
#
#import glob
#shearkey="shear-only-BIAS"
#files=glob.glob("*%s*0.*.txt"%shearkey)
#files.sort()
#post_shear = []
#x=[]
#bias_sig8_s = []
#biasx=[]
#if shearkey:
#    for i, f in enumerate(files):
#        post_shear+=[np.loadtxt(f).T[1]]
#        x+=[np.loadtxt(f).T[0]]
#        m = np.argwhere( np.exp(post_shear[-1])==np.exp(post_shear[-1]).max() )[0]
#        bias_sig8_s += [ x[-1][m][0] -0.82 ]
#        v = [float(t) for t in f.split('-') if ("." in t and "0" in t)][0]
#        biasx += [v]
#
#sgkey="shear+ggl-BIAS"
#files=glob.glob("*%s*0.*.txt"%sgkey)
#files.sort()
#post_sg = []
#xsg=[]
#bias_sig8_sg = []
#biasx_sg=[]
#if sgkey:
#    for i, f in enumerate(files):
#        post_sg+=[np.loadtxt(f).T[1]]
#        xsg+=[np.loadtxt(f).T[0]]
#        m = np.argwhere( np.exp(post_sg[-1])==np.exp(post_sg[-1]).max() )[0]
#        bias_sig8_sg += [ xsg[-1][m][0] -0.82 ]
#        v = [float(t) for t in f.split('-') if ("." in t and "0" in t)][0]
#        biasx_sg += [v]
#
#sgpkey="shear+ggl+pos-BIAS"
#files=glob.glob("*%s*0.*.txt"%sgpkey)
#files.sort()
#post_sgp = []
#xsgp=[]
#bias_sig8_sgp = []
#biasx_sgp=[]
#if sgpkey:
#    for i, f in enumerate(files):
#        post_sgp+=[np.loadtxt(f).T[1]]
#        xsgp+=[np.loadtxt(f).T[0]]
#        m = np.argwhere( np.exp(post_sgp[-1])==np.exp(post_sgp[-1]).max() )[0]
#        bias_sig8_sgp += [ xsgp[-1][m][0] -0.82 ]
#        v = [float(t) for t in f.split('-') if ("." in t and "0" in t)][0]
#        biasx_sgp += [v]
#
#spkey="shear+pos-BIAS"
#files=glob.glob("*%s*0.*.txt"%spkey)
#files.sort()
#post_sp = []
#xsp=[]
#bias_sig8_sp = []
#biasx_sp=[]
#if spkey:
#    for i, f in enumerate(files):
#        post_sg+=[np.loadtxt(f).T[1]]
#        xsp+=[np.loadtxt(f).T[0]]
#        m = np.argwhere( np.exp(post_sp[-1])==np.exp(post_sp[-1]).max() )[0]
#        bias_sig8_sp += [ xsp[-1][m][0] -0.82 ]
#        v = [float(t) for t in f.split('-') if ("." in t and "0" in t)][0]
#        biasx_sp += [v]
#
#plt.plot(biasx, np.array(bias_sig8_s)/0.82, 'm-', label='WL')
#plt.plot(biasx, np.array(bias_sig8_s)/0.82, 'mo')
#plt.plot(biasx_sg, np.array(bias_sig8_sg)/0.82, 'g-', label='WL+ggl')
#plt.plot(biasx_sg, np.array(bias_sig8_sg)/0.82, 'go')
#plt.plot(biasx_sgp, np.array(bias_sig8_sgp)/0.82, 'b-', label='WL+ggl+LSS')
#plt.plot(biasx_sgp, np.array(bias_sig8_sgp)/0.82, 'bo')
#
#plt.xlabel('bias $\delta z$')
#plt.ylabel('bias $\Delta \sigma_8 / \sigma_8$')
#plt.legend(loc='upper right')
