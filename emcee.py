import numpy as np
import scipy as sp
import astropy.io.fits as pyfits
import astropy.table as tb
import os, pdb, copy
import pylab as plt
from tools import samplers as samp

labels={"sigma8": "$\sigma_8$", "s8" : "$S_8 \equiv \sigma_8 ( \Omega_m / 0.31 ) ^{0.5}$", "w0":"$w_0$", "omega_m": "$\Omega _m$", "h" : '$h$', "ns":"$n_s$", "omega_b" : '$\Omega _b$', "shear--bias_1":"$\delta z^1$", "shear--bias_2":"$\delta z^2$", "shear--bias_3":"$\delta z^3$", "sz1":"$S_z ^1$", "sz2":"$S_z ^2$", "sz3":"$S_z ^3$", "redmagic--bias_1":"$\delta z^1$ (redmagic)", "redmagic--bias_2":"$\delta z^2$ (redmagic)", "redmagic--bias_3":"$\delta z^3$ (redmagic)", "a": "$A$", "eta": "$\eta$","A_II": "$A_{II}$", "eta_II": "$\eta_{II}$", "A_GI": "$A_{GI}$", "eta_GI": "$\eta_{GI}$", "m1":"$m_1$", "m2":"$m_2$", "m3":"$m_3$", "b_1": '$b_g^1$', "b_2" : '$b_g^2$', "b_3" : '$b_g^3$'}


parameters={"omega_m" : 'cosmological_parameters--omega_m', "s8":"cosmological_parameters--s8", "w0": "cosmological_parameters--w", "h" : 'cosmological_parameters--h0', "omega_b" : 'cosmological_parameters--omega_b', "sigma8": 'cosmological_parameters--sigma8_input', "ns": 'cosmological_parameters--n_s', "A_II" : 'intrinsic_alignment_parameters--a_ii', "eta_II" : 'intrinsic_alignment_parameters--alpha_ii', "A_GI" : 'intrinsic_alignment_parameters--a_gi', "eta_GI" : 'intrinsic_alignment_parameters--alpha_gi', "A" : 'intrinsic_alignment_parameters--a',  "eta" : 'intrinsic_alignment_parameters--alpha', "b_1": 'bias_parameters--b_1', "b_2" : 'bias_parameters--b_2', "b_3" : 'bias_parameters--b_3', "shear--bias_1": 'shear--bias_1', "shear--bias_2" : 'shear--bias_2', "shear--bias_3" : 'shear--bias_3',  "sz1": 'shear--s_z_1', "sz2" : 'shear--s_z_2', "sz3" : 'shear--s_z_3', "m1" : 'shear--m1', "m2" : 'shear--m2', "m3" : 'shear--m3', "redmagic--bias_1" : 'redmagic--bias_1', "redmagic--bias_2" : 'redmagic--bias_2', "redmagic--bias_3" : 'redmagic--bias_3'}

fiducial = {"sigma8":0.82, "s8":0.8273733018687 ,"w0": -1.0, "omega_m":0.3156, "eta":0.0 ,"eta_ii": 0.0, "eta_gi": 0.0, "a_ii": 1.0 , "a_gi": 1.0, "a":1.0}

class chain(samp.sampler):
	def __init__(self, filename):
		self.samples=tb.Table.read(filename, format="ascii")
		try:
			self.post = self.samples["post"]
			self.samples.remove_column("post")
		except:
			self.post = self.samples["like"]
			self.samples.remove_column("like")
		self.wt = np.ones_like(self.post)
		

		self.mask = np.ones_like(self.post)
		self.bounds={}

		self.filename=filename

		print "found %d posterior samples in %s"%(len(self.post), self.filename)

		sep="END_OF_PRIORS_INI\n"
		text = open(self.filename).read()
		self.header = text.split(sep)[0]+sep
		self.npar = int(self.header.split("n_varied=")[1].split("\n")[0])

	def add_column(self, name, pmin=-1, pmax=1, values="random", cosmosis_section="cosmological_parameters"):
		"""Draw new samples from a uniform distribution over a specified range."""
		
		num = len(self.samples)
		if values=="random":
			newcol = (np.random.rand(num)*(pmax-pmin)) + (pmin+pmax)/2. - (pmax-pmin)/2
			print "Added random samples "
			print newcol
		# A rather hacked together bit of code to parse a human readable string
		# specifying a specific combination of parameters	
		else:
			symbols = ["*","/","+","-","**","(",")"]
			par = values
			for s in symbols:
				par=par.replace(s, " ")

			par_names = []
			for p in par.split(" "):
				try:
					v=float(p)
				except:
					if not p: 
						continue
					else:
						par_names.append(p)

			for p in par_names:
				values = values.replace(p, "self.samples[parameters['%s']]"%p)

			print "Generating new column %s from existing columns "%name , par_names
			print values 
			exec "newcol = " + values

		newcol = tb.Column(newcol, name=cosmosis_section+"--"+name)


		self.samples.add_column(newcol, index=len(self.samples.dtype))

		self.header = self.header.replace("\tpost", "\t%s--%s\tpost"%(cosmosis_section, name))

		self.header = self.header.replace("n_varied=%d"%self.npar, "n_varied=%d"%(self.npar+1))
		self.npar+=1

	def add_weights(self, weights):
		"""Store an array of weights."""
		self.wt=weights

	def burn(self, num):
		"""Discard a specified number of samples from the start of the chain."""
		
		n0 = len(self.samples)

		self.samples = self.samples[num:]
		self.post = self.post[num:]

		self.remove_bounds()

		print "Discarding %d/%d samples"%(num, n0)

	def thin(self, factor):
		print "Thinning chain by a factor of %2.3f"%factor
		num  = len(self.samples[range[0]:range[1]])/factor
		sel = np.arange(0,len(self.samples), len(self.samples)/num).astype(int)

		self.samples = self.samples[sel]
		self.post = self.post[sel]

	def change_deltaz_prior(self, newprior, oldprior):
		b1 = self.samples["shear--bias_1"]
		b2 = self.samples["shear--bias_2"]
		b3 = self.samples["shear--bias_3"]

		pold = gaussian(b1,0.0,oldprior)*gaussian(b2,0.0,oldprior)*gaussian(b3,0.0,oldprior)
		pnew = gaussian(b1,0.0,newprior)*gaussian(b2,0.0,newprior)*gaussian(b3,0.0,newprior)

		self.wt = pnew/pold



	def do_importance_resampling(self):
		newsamp=[]
		for i,samp in enumerate(self.samples):
			draw=np.random.rand()
			if draw<self.wt[i]:
				newsamp.append(i)

		nsamp=len(newsamp)

		print "New chain contains %d/%d points after resampling"%(len(newsamp), len(self.samples))

		return nsamp, np.array(newsamp), self.samples[np.array(newsamp)]

	def write_columns(self, filename=None, overwrite=False, bounds=True, impweight=False, apply_weights=False, halfchains=False):
		if not filename and overwrite:
			filename = self.filename
			print "Overwriting file %s"%self.filename

		if bounds:
			sel = np.array(self.mask).astype(bool)
		else:
			sel = np.ones_like(self.post).astype(bool)

		samp = copy.deepcopy(self.samples)
		nsamp=len(self.samples)

		if impweight:
			print "Using importance sampled posterior ",np.array(self.wt)
			nsamp, i, samp = self.do_importance_resampling()
			print i
			post=np.array(self.post)[i]
			sel=np.ones_like(post).astype(bool)
		elif apply_weights:
			print "Applying weights to posterior ",self.wt
			post=self.post*self.wt
		else:
			post=self.post
		post = tb.Column(post, name="post")
		samp.add_column(post)

		np.savetxt(filename, np.array(samp[sel]).T, header=self.header, comments="")

		if halfchains:
			import os
			base=os.path.basename(filename)
			dirname=os.path.dirname(filename)
			h1=os.path.join(dirname,"h1_"+base)
			h2=os.path.join(dirname,"h2_"+base)
			self.write_halfchains(0,name1=h1,name2=h2)

		del(samp)

		return nsamp

		#self.samples.remove_column("post")

	def write_halfchains(self, burn, name1="halfchain1.txt", name2="halfchain2.txt"):
		post = tb.Column(self.post, name="post")
		samp = copy.deepcopy(self.samples)
		samp.add_column(post)

		n = len(self.samples)

		np.savetxt(name1, np.array(samp[:n/2]).T, header=self.header, comments="")
		np.savetxt(name2, np.array(samp[n/2:]).T, header=self.header, comments="")

		print "Saving half chains to %s %s"%(name1, name2)

		del(samp)

	def show_convergence(self, param, colour="m.", range=(0,-1), thin=None, offset=0, newplot=False):
		if newplot:
			plt.figure(newplot)
		if thin:
			print "Thinning chain by a factor of %2.3f"%thin
			num  = len(self.samples[range[0]:range[1]])/thin
			sel = np.arange(0,len(self.samples[range[0]:range[1]]), len(self.samples)/num).astype(int)
		else:
			sel = np.ones_like(self.samples[range[0]:range[1]]).astype(bool)
		plt.plot(self.samples[parameters[param]][range[0]:range[1]][sel]+offset, colour)

	def add_rows(self, num):
		newrows = np.zeros((len(self.samples.dtype.names), num))
		newrows = tb.Table(newrows.T, names=self.samples.dtype.names)

		self.samples = tb.vstack((newrows, self.samples))

		self.post = np.hstack((np.zeros(num), self.post))

		print "Added %d dummy samples."%num

	def check_convergence(self, param, burn):
		samples = self.samples[parameters[param]][burn:]
		post = self.post[burn:]

		print "Mean of samples: %3.4f"%np.mean(samples)
		median = np.median(samples)
		print "Median of samples: %3.4f"%median
		best_fit = samples[np.argwhere( post==post.max() )[0,0]]
		print "Max. likelihood value: %3.4f"%best_fit
		print "Standard deviation of samples: %3.4f"%np.std(samples)

	def smooth(self, plot_type="variance",param="s8", nbins=50, line=True, colour="m"):
		samp=self.samples[parameters[param]]
		bins=np.linspace(0,len(samp),nbins).astype(int)
		print bins
		x=(bins[1:]+bins[:-1])/2.

		y=[]

		for i,lower in enumerate(bins[:-1]):
			upper=bins[i+1]
			if plot_type=="variance": y.append(np.std(samp[lower:upper]))
			elif plot_type=="mean": y.append(np.mean(samp[lower:upper]))

		plt.plot(x,y,".", color=colour)
		if line:
			plt.plot(x,y,"-", color=colour)



	def plot_statistic_vs_burnin(self, param, points=20, mean=True, error=False, showpoints=True):
		samples = self.samples[parameters[param]]
		burn = np.linspace(0, len(samples), points)[:-1].astype(int)
		means = []
		stds=[]
		for b in burn:
			if mean:
				means.append( np.mean(samples[b:]) )

			if error:
				stds.append(np.std(samples[b:]))


		if mean:
			plt.plot(burn, means, "-", color="purple", lw=2.0, label="mean")
			if showpoints:
				plt.plot(burn, means, "o", color="purple")
			plt.ylabel("mean %s"%labels[param])

		if error:
			plt.plot(burn, stds, "-", color="midnightblue", lw=2.0, label="error")
			if showpoints:
				plt.plot(burn, stds, "o", color="midnightblue")
			plt.ylabel("error %s"%labels[param])

		if error and mean:
			plt.ylabel("error or mean %s"%labels[param])
			plt.legend(loc="upper left")


		plt.xlabel("points discarded")
		

		plt.axvline(len(samples), color="k", lw=2)


	def hist(self, param, bins=60, mark_fiducial=True, w=None, bounds=False, mean=False):
		p1 = parameters[param]
		plt.xlabel(labels[param])
		if bounds:
			samp = self.samples[p1][np.array(self.mask.astype(bool))]
		else:
			samp = self.samples[p1]

		plt.hist(samp, bins=bins, weights=w, histtype="step", normed=True)
		if mark_fiducial:
			plt.axvline(fiducial[param], color="k", ls="-.")
		if mean:
			plt.axvline(np.mean(samp), color="k", ls=":")

	def set_bounds(self, param, limits):
		self.bounds[param] = limits
		sel = (np.array(self.samples[parameters[param]])>limits[0]) & (np.array(self.samples[parameters[param]])<limits[1])
		self.mask *= sel.astype(int)
		print "Imposing bounds %f < %s < %f"%(limits[0], param, limits[1])

	def remove_bounds(self):
		self.bounds = {}
		self.mask=np.ones_like(self.post)

	def scatter(self, param=["sigma8", "omega_m"], colour=None, ymin=None, ymax=None, xmin=None, xmax=None, burn=0, bounds=True, mark_fiducial=False):
		p1 = parameters[param[0]]
		p2 = parameters[param[1]]

		sel1 = np.linspace(0, len(self.samples[p1]),len(self.samples[p1]))  > burn

		#Impose parameter bounds if required
		if bounds:
			sel2 = self.mask.astype(bool)
		else:
			sel2 = np.ones_like(self.post).astype(bool)
		print len(self.samples[p1])
		print len(self.samples[p1][sel1 & sel2])
		samp1 = self.samples[p1][sel1 & sel2]
		samp2 = self.samples[p2][sel1 & sel2]

		if not colour:
			colourbar=True
			colour = self.post[sel1 & sel2]

		else:
			colour=False

		plt.scatter(samp1, samp2, marker=".", c=colour)

		plt.xlabel(labels[param[0]])
		plt.ylabel(labels[param[1]])

		if mark_fiducial:
			print fiducial[param[0]], fiducial[param[1]]
			plt.plot(fiducial[param[0]], fiducial[param[1]], "x", color='k', markersize=10, markerfacecolor='k')

		if colourbar:
			plt.colorbar(label="posterior")

	def check_compatibility(self, chain):
		compatible = True
		if chain.npar!=self.npar:
			print "parameter space dimensons differ"
			compatible = False

		if not np.in1d(self.samples.dtype.names, chain.samples.dtype.names).all():
			print "parameter space is non identical"
			compatible - False

		return compatible

def merge_chains(chain1, chain2):
	if not chain1.check_compatibility(chain2):
		raise ValueError("Chains are not compatible.")
	else:
		new = copy.deepcopy(chain1)
		new.samples = np.hstack((chain1.samples, chain2.samples))
		new.post = np.hstack((chain1.post, chain2.post))

	return new

def gaussian(x,mu,sigma):
	return np.exp(-1.0*(x-mu)*(x-mu)/2/sigma/sigma) / np.sqrt(2*np.pi)/sigma


def extract_confidence_bounds(directory, param="s8", conf="std", deltaz=False):
	if isinstance(conf, str):
		if conf.lower()=="std":
			print "Using standard deviation of points as the errorbar on %s"%param
			means=open(directory+"/means.txt").read().split("\n")
			y=[]
			dz=[]
			for line in means:
				if deltaz:
					if ".txt" in line:
						print line.split("_deltazprior")
						if "_deltazprior" in line:
							val=line.split("_deltazprior")[1].split("_")[0]
						elif "shear_" not in line:
							val=line.split("_")[1].replace(".txt","")
						else:
							val=line.split("shear_")[1].replace(".txt","")

						#except:
							#val = line.split("_")[-1].replace(".txt", "").replace("postprocessed","")
						dz.append(float(val))
				try:
					if param in line.split()[0]:
						print line.split()
						y.append(float(line.split()[2]))

				except:
					continue

			if not deltaz:
				return np.array(y)
			else:
				return dz, np.array(y)






	lower=open(directory+"/low%d.txt"%conf).read().split("\n")
	upper=open(directory+"/upper%d.txt"%conf).read().split("\n")

	l=[]
	dz=[]
	for line in lower:
		try:
			if param in line.split()[0]:
				l.append(float(line.split()[1]))

		except:
			continue

	u=[]
	for line in upper:
		print line
		if deltaz:
			if "deltazprior" in line:
				val = line.split("_deltazprior")[1].split("_")[0]
				dz.append(float(val))
		try:
			if param in line.split()[0]:
				u.append(float(line.split()[1]))
		except:
			continue

	print "Upper bounds", u
	print "Lower bounds", l

	if not deltaz:
		return np.array(u)-np.array(l)
	else:
		return dz, np.array(u)-np.array(l)


def extract_bias(directory, param="s8", deltaz=False):
	bias=open(directory+"/means.txt").read().split("\n")

	dz=[]
	b=[]
	for line in bias:
		print line
		if deltaz:
			if "deltazprior" in line:
				dz.append(float(line.replace(".txt", "").split("_deltazprior")[1].split("_")[0].split("/")[0]))
		try:
			if param in line.split()[0]:
				b.append(float(line.split()[1]))
		except:
			continue

	print "Means:", b

	b=np.array(b)-fiducial[param]

	if not deltaz:
		return b
	else:
		return dz, b



# -*- coding: utf-8 -*-


__all__ = ["function", "integrated_time", "AutocorrError"]


def function(x, axis=0, fast=False):
    """Estimate the autocorrelation function of a time series using the FFT.
    Args:
        x: The time series. If multidimensional, set the time axis using the
            ``axis`` keyword argument and the function will be computed for
            every other axis.
        axis (Optional[int]): The time axis of ``x``. Assumed to be the first
            axis if not specified.
        fast (Optional[bool]): If ``True``, only use the first ``2^n`` (for
            the largest power) entries for efficiency. (default: False)
    Returns:
        array: The autocorrelation function of the time series.
    """
    x = np.atleast_1d(x)
    m = [slice(None), ] * len(x.shape)

    # For computational efficiency, crop the chain to the largest power of
    # two if requested.
    if fast:
        n = int(2**np.floor(np.log2(x.shape[axis])))
        m[axis] = slice(0, n)
        x = x
    else:
        n = x.shape[axis]

    # Compute the FFT and then (from that) the auto-correlation function.
    f = np.fft.fft(x - np.mean(x, axis=axis), n=2*n, axis=axis)
    m[axis] = slice(0, n)
    acf = np.fft.ifft(f * np.conjugate(f), axis=axis)[m].real
    m[axis] = 0
    return acf / acf[m]


def integrated_time(x, low=10, high=None, step=1, c=10, full_output=False,
                    axis=0, fast=False):
    """Estimate the integrated autocorrelation time of a time series.
    This estimate uses the iterative procedure described on page 16 of `Sokal's
    notes <http://www.stat.unc.edu/faculty/cji/Sokal.pdf>`_ to determine a
    reasonable window size.
    Args:
        x: The time series. If multidimensional, set the time axis using the
            ``axis`` keyword argument and the function will be computed for
            every other axis.
        low (Optional[int]): The minimum window size to test. (default: ``10``)
        high (Optional[int]): The maximum window size to test. (default:
            ``x.shape[axis] / (2*c)``)
        step (Optional[int]): The step size for the window search. (default:
            ``1``)
        c (Optional[float]): The minimum number of autocorrelation times
            needed to trust the estimate. (default: ``10``)
        full_output (Optional[bool]): Return the final window size as well as
            the autocorrelation time. (default: ``False``)
        axis (Optional[int]): The time axis of ``x``. Assumed to be the first
            axis if not specified.
        fast (Optional[bool]): If ``True``, only use the first ``2^n`` (for
            the largest power) entries for efficiency. (default: False)
    Returns:
        float or array: An estimate of the integrated autocorrelation time of
            the time series ``x`` computed along the axis ``axis``.
        Optional[int]: The final window size that was used. Only returned if
            ``full_output`` is ``True``.
    Raises
        AutocorrError: If the autocorrelation time can't be reliably estimated
            from the chain. This normally means that the chain is too short.
    """
    size = 0.5 * x.shape[axis]
    if int(c * low) >= size:
        raise AutocorrError("The chain is too short")

    # Compute the autocorrelation function.
    f = function(x, axis=axis, fast=fast)

    # Check the dimensions of the array.
    oned = len(f.shape) == 1
    m = [slice(None), ] * len(f.shape)

    # Loop over proposed window sizes until convergence is reached.
    if high is None:
        high = int(size / c)
    for M in np.arange(low, high, step).astype(int):
        # Compute the autocorrelation time with the given window.
        if oned:
            # Special case 1D for simplicity.
            tau = 1 + 2 * np.sum(f[1:M])
        else:
            # N-dimensional case.
            m[axis] = slice(1, M)
            tau = 1 + 2 * np.sum(f[m], axis=axis)

        # Accept the window size if it satisfies the convergence criterion.
        if np.all(tau > 1.0) and M > c * tau.max():
            if full_output:
                return tau, M
            return tau

        # If the autocorrelation time is too long to be estimated reliably
        # from the chain, it should fail.
        if c * tau.max() >= size:
            break

    raise AutocorrError("The chain is too short to reliably estimate "
                        "the autocorrelation time")


class AutocorrError(Exception):
    """Raised if the chain is too short to estimate an autocorrelation time.
    """
    pass
