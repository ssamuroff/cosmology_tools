import numpy as np
import scipy as sp
import astropy.io.fits as pyfits
import os, pdb
from tools import fisher as fi
from samplers import sampler
import tools.emcee as mc
import tools.diagnostics as di
import tools.arrays as arr
import fitsio

import matplotlib.colors
import matplotlib
import pylab as plt
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['font.size']=14
matplotlib.rcParams['legend.fontsize']=14
matplotlib.rcParams['xtick.major.size'] = 10.0
matplotlib.rcParams['ytick.major.size'] = 10.0

colours=['midnightblue','forestgreen','pink', 'purple', "lightblue"]
col = {"wl" : "midnightblue", "wl+lss" : 'forestgreen', "ggl+lss" : 'pink',  "wl+ggl+lss" : 'purple',  "wl+ggl" : "lightblue"}


param_names={"sigma8_input" : ["cosmological_parameters--sigma8_input", "\sigma_8"], "s8":["cosmological_parameters--s8", "S_8"], "omega_m": ["cosmological_parameters--omega_m", "\Omega_m"],  "w" : ["cosmological_parameters--w", "w_0"], "wa" : ["cosmological_parameters--wa", "w_a"]}

fiducial={"sigma8_input" : 0.82, "omega_m":0.312, "w" : -1.0, "wa" : 0.0, "s8":0.8273733018687}

probe_labels={"wl": "\gamma \gamma", "ggl":"\gamma \delta _g", "lss" : "\delta _g \delta _g", "wl+ggl+lss":"\gamma \gamma + \gamma \delta _g + \delta _g \delta _g" }

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



	def make(self, output="show", loc=None, normalise=False, linestyle="-", xlim_upper=0.06, ylim_lower=None, ylim_upper=None, s=True, sp=True, gp=True, sg=False, sgp=True, param="sigma8_input", extra_label="", colour=None):
		import pylab as plt
		import matplotlib.colors
		import matplotlib
		matplotlib.rcParams['font.family']='serif'
		matplotlib.rcParams['font.size']=16

		if isinstance(param, list):
			name="nd_contour_area"
			for p in param:
				name+="--%s"%p
		else:
			name = param_names[param][0]

		if extra_label:
			extra_label=" "+extra_label

		if colour is None:
			col0=col["wl+ggl+lss"]
		else:
			col0=colour

		if normalise:
			if s:
				plt.plot(self.shear_shear[0], self.shear_shear[1][name]/self.shear_shear[1][name][0], col["wl"], label="$%s$"%probe_labels["wl"]+extra_label, linewidth=2.0, linestyle=linestyle)
			if sp:
				plt.plot(self.sp[0], self.sp[1][name]/self.sp[1][name][0], col["wl+lss"], label="$%s+%s$"%(probe_labels["wl"],probe_labels["lss"])+extra_label, linewidth=2.0, linestyle=linestyle)
			if gp:
				plt.plot(self.gp[0], self.gp[1][name]/self.gp[1][name][0], col["ggl+lss"], label="%s+%s"%(probe_labels["ggl"],probe_labels["lss"])+extra_label, linewidth=2.0, linestyle=linestyle)
			if sg:
				plt.plot(self.sg[0], self.sg[1][name]/self.sg[1][name][0], col["wl+ggl"], label="$%s+%s$"%(probe_labels["wl"],probe_labels["ggl"])+extra_label, linewidth=2.0, linestyle=linestyle)
			if sgp:
				plt.plot(self.sgp[0], self.sgp[1][name]/self.sgp[1][name][0], col0, label="$%s+%s+%s$"%(probe_labels["wl"],probe_labels["ggl"],probe_labels["lss"])+extra_label, linewidth=2.0, linestyle=linestyle)
			plt.ylabel("error degradation $\Delta %s/ \Delta %s (\Delta \delta z=0)$"%(param_names[param][1], param_names[param][1]), fontsize=20)

		else:
			norm = abs(fiducial[param])
			if norm==0.0:
				norm=1.0
				ab=True
			else:
				ab=False
			if s:
				plt.plot(self.shear_shear[0], self.shear_shear[1][name]/norm, col["wl"], label="$%s$"%probe_labels["wl"]+extra_label, linewidth=2.0, linestyle=linestyle)
			if sp:
				plt.plot(self.sp[0], self.sp[1][name]/norm, col["wl+lss"], label="$%s+%s$"%(probe_labels["wl"],probe_labels["lss"])+extra_label, linewidth=2.0, linestyle=linestyle)
			if gp:
				plt.plot(self.gp[0], self.gp[1][name]/norm, col["ggl+lss"], label="$%s+%s$"%(probe_labels["ggl"],probe_labels["lss"])+extra_label, linewidth=2.0, linestyle=linestyle)
			if sg:
				plt.plot(self.sg[0], self.sg[1][name]/norm, col["wl+ggl"], label="$%s+%s$"%(probe_labels["wl"],probe_labels["ggl"])+extra_label, linewidth=2.0, linestyle=linestyle)
			if sgp:
				plt.plot(self.sgp[0], self.sgp[1][name]/norm, col0, label="$%s+%s+%s$"%(probe_labels["wl"],probe_labels["ggl"], probe_labels["lss"])+extra_label, linewidth=2.0, linestyle=linestyle)
			if not ab: 
				plt.ylabel("fractional error $\Delta %s / %s$"%(param_names[param][1], param_names[param][1]), fontsize=20)
			else:
				plt.ylabel("absolute error $\Delta %s$"%(param_names[param][1]), fontsize=20)


		plt.xlabel("prior width $\Delta \delta z$", fontsize=20)
		plt.legend(loc="lower right")
		plt.xlim(xmax=xlim_upper)
		plt.ylim(ymin=ylim_lower, ymax=ylim_upper)
		plt.tight_layout()

		if output=="show":
			plt.show()
		elif output=="save":
			plt.savefig(loc)
		elif output=="none":
			return 0


		plt.close()



	def make2d(self, loc=None, s=True, gp=False, sp=True, sgp=True, normalise=False, xlim_upper=0.06, ylim_lower=None, ylim_upper=None, param=["sigma8_input", "omega_m"]):
		import pylab as plt

		p1 = param[0]
		p2 = param[1]

		if normalise:
			if s:
				plt.plot(self.shear_shear_2d[0], self.shear_shear_2d[1][(p1,p2)]/self.shear_shear_2d[1][(p1,p2)][0], colours[0], label="WL", linewidth=2.0)
			if sp:
				plt.plot(self.sp_2d[0], self.sp_2d[1][(p1,p2)]/self.sp_2d[1][(p1,p2)][0], colours[1], label="WL+LSS", linewidth=2.0)
			if gp:
				plt.plot(self.gp_2d[0], self.gp_2d[1][(p1,p2)]/self.gp_2d[1][(p1,p2)][0], colours[2], label="GGL+LSS", linewidth=2.0)
			if sgp:
				plt.plot(self.sgp_2d[0], self.sgp_2d[1][(p1,p2)]/self.sgp_2d[1][(p1,p2)][0], colours[3], label="WL+GGL+LSS", linewidth=2.0)
			#pdb.set_trace()
			plt.ylabel("error degradation $FOM^{-1}_{%s %s} / FOM^{-1}_{%s %s} (\Delta \delta z=0)$"%(param_names[param[0]][1], param_names[param[1]][1], param_names[param[0]][1], param_names[param[1]][1]))
		else:
			if s:
				plt.plot(self.shear_shear_2d[0], self.shear_shear_2d[1][(p1,p2)], colours[0], label="WL", linewidth=2.0)
			if sp:
				plt.plot(self.sp_2d[0], self.sp_2d[1][(p1,p2)], colours[1], label="WL+LSS", linewidth=2.0)
			if gp:
				plt.plot(self.gp_2d[0], self.gp_2d[1][(p1,p2)], colours[2], label="GGL+LSS", linewidth=2.0)
			if sgp:
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




class fig9:
	def __init__(self, postprocess_dir=".", sec=None, sp=False, sg=False, sgp=True, gp=False, from_fisher=False, prior=None, param="sigma8_input"):
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
			self.process_fisher(from_fisher, prior, sp, sgp, gp, sg, sec, param )

	def process_fisher(self, filenames, prior, sp=False, sgp=True, gp=False, sg=False, sec=None, param="sigma8_input"):
		self.shear_shear=[]
		self.sgp=[]
		self.sp=[]
		self.gp=[]
		self.sg=[]
		for filename in filenames:
			f = fi.fisher(filename)
			
			if not sec:
				sec = "nz_shear_errors--bias"

			nparam = str(list(np.unique(f.param))).count(sec)
			sections=[]
			for i in xrange(nparam):
				sections.append(sec+"_%d"%(i+1))
			f.remove_priors(sections, np.array(nparam*[prior]))

			# If sigma8 isn't a free parameter try omegam instead
			if param not in f.param:
				param="omega_m"
			x, y = fi.photoz_prior_width_plot(f, bias_name=sec, parameter=param)

			ngal = float(filename.split("_")[2].replace("ngal", ""))
			print "ngal:%f"%ngal

			err = {param_names[param][0]:y}
			

			if "shear+ggl+pos" in filename:
				self.sgp += [[ngal, x, err]]
			if "shear+pos" in filename:
				self.sp += [[ngal, x, err]]
			if "ggl+pos" in filename:
				self.gp += [[ngal, x, err]]
			if ("shear+ggl" in filename) and ("pos" not in filename):
				self.sg += [[ngal, x, err]]
			if ("ggl" not in filename) and ("pos" not in filename):
				self.shear_shear += [[ngal, x, err]]

	def reshape_to_grid(self, sp=False, sgp=True, gp=False, sg=False, param="sigma8_input"):
		self.coarse_grid = {}
		if sgp:
			self.sgp = np.array(self.sgp)
			order=np.argsort(self.sgp.T[0])
			self.sgp.T[0]=self.sgp.T[0][order]
			self.sgp.T[1]=self.sgp.T[1][order]
			self.sgp.T[2]=self.sgp.T[2][order]

			self.grid_sgp = []
			for i, sgp in enumerate(self.sgp):
				self.grid_sgp.append(sgp[2][param_names[param][0]])

			self.grid_sgp=np.array(self.grid_sgp)

		if sp:
			self.sp = np.array(self.sp)
			order=np.argsort(self.sp.T[0])
			self.sp.T[0]=self.sp.T[0][order]
			self.sp.T[1]=self.sp.T[1][order]
			self.sp.T[2]=self.sp.T[2][order]
			self.grid_sp = []
			for i, sp in enumerate(self.sp):
				self.grid_sp.append(sp[2][param_names[param][0]])
			self.grid_sp=np.array(self.grid_sp)

		if gp:
			order=np.argsort(self.gp.T[0])
			self.gp.T[0]=self.gp.T[0][order]
			self.gp.T[1]=self.gp.T[1][order]
			self.gp.T[2]=self.gp.T[2][order]
			self.gp = np.array(self.gp)
			self.grid_gp = []

			for i, gp in enumerate(self.gp):
				self.grid_gp.append(gp[2][param_names[param][0]])

			self.grid_gp=np.array(self.grid_gp)

		if sg:
			order=np.argsort(self.sg.T[0])
			self.sg.T[0]=self.sg.T[0][order]
			self.sg.T[1]=self.sg.T[1][order]
			self.sg.T[2]=self.sg.T[2][order]
			self.sg = np.array(self.gp)
			self.grid_sg = []

			for i, sg in enumerate(self.sg):
				self.grid_sg.append(sg[2][param_names[param][0]])

			self.grid_sg=np.array(self.grid_sg)

	def make(self, sp=False, sgp=True, gp=False, sg=False, param="sigma8_input", title="WL+GGL+LSS", normalise=False):

		self.reshape_to_grid(sp=sp, sgp=sgp, gp=gp, sg=sg, param=param)

		import scipy.interpolate as inter
		import pylab as plt
		ngal = self.sgp.T[0]
		if len(ngal)>20:
			interp_kind="cubic"
		else:
			interp_kind="linear"
		b = self.sgp.T[1][0]

		print "Using %s interpolation"%interp_kind

		f = inter.interp2d(b, ngal, self.grid_sgp, kind=interp_kind)
		xf = np.linspace(1e-50, 0.06,200)
		yf = np.linspace(0.01,4.5,200)

		grid = f(xf,yf)

		if not normalise:
			norm = fiducial[param]
		else:
			norm = grid[-1,0]

		plt.imshow(grid/norm, extent=[xf[0],xf[-1],yf[0],yf[-1]], aspect=0.01, origin="lower")
		if not normalise:
			plt.colorbar(label="Statistical error $\Delta \sigma_8 / \sigma_8$")
		else:
			plt.colorbar(label="Error degradation $\Delta \sigma_8 / \Delta \sigma_{8,0}$")

		plt.title(title)

		CS = plt.contour(xf, yf, grid/norm, 5, colors="k")
		plt.clabel(CS,inline=1, fontsize=10, fmt="%3.6f")
		plt.xlabel("prior width $\Delta \delta z_{LSS}$")
		plt.ylabel("number density of LSS sample (arcmin$^{-2}$)")
		plt.show()

class fig1:
	def __init__(self, filename):
		import cPickle as pi
		print "Getting data from :%s"%filename
		self.data = pi.load((open(filename)))
	def make(self, nbin_wl=3, nbin_lss=3):
		import matplotlib.colors
		import matplotlib
		matplotlib.rcParams['font.family']='serif'
		matplotlib.rcParams['font.size']=16
		matplotlib.rcParams['legend.fontsize']=15
		matplotlib.rcParams['xtick.major.size'] = 10.0
		matplotlib.rcParams['ytick.major.size'] = 10.0

		plt.subplots_adjust(wspace=0, hspace=0)
		xs = self.data["ell_shear"]
		xg = self.data["ell_ggl"]
		xp = self.data["ell_pos"]
		for j in xrange(1,nbin_lss+1):
			for i in xrange(1,nbin_wl+1):
				plot_index = i + (j-1)*3
				print i, j, plot_index
				ax=plt.subplot(nbin_lss, nbin_wl, plot_index )
				try:
					plt.plot(xs, xs*xs*self.data["C_ee"]["bin_%d_%d"%(i,j)], "darkviolet", linestyle="-", lw=2.0)
				except:
					    print "No shear Cl"
				try:
					plt.plot(xg, xg*xg*self.data["C_ne"]["bin_%d_%d"%(i,j)], "forestgreen", linestyle="-.", lw=2.0)
				except:
					print "No ggl Cl"
				try:
					plt.plot(xp, xp*xp*self.data["C_nn"]["bin_%d_%d"%(i,j)], "darkred", linestyle=":", lw=2.0)
				except:
					print "No position Cl"
				plt.xscale("log")
				plt.yscale("log")
				plt.xlim(10., 2700)
				plt.ylim(1.0e-7, 3.0)
				plt.annotate("%d,%d"%(i,j), xy=(1000.0, 6.0e-7), fontsize=15)
				if i!=1:
					plt.setp(ax.get_yticklabels(), visible=False)
				elif plot_index==1:
					ax.set_yticks([1.0e-6, 1.0e-4, 1.0e-2, 1.0])
					ax.set_yticklabels(["$10^{-6}$", "$10^{-4}$", "$10^{-2}$", "$10^{0}$"], fontsize=18)
				else:
					ax.set_yticks([1.0e-6, 1.0e-4, 1.0e-2, 1.0])
					ax.set_yticklabels(["$10^{-6}$", "$10^{-4}$", "$10^{-2}$", "$10^{0}$"], fontsize=18)
				if j!=nbin_lss:
					plt.setp(ax.get_xticklabels(), visible=False)
				else:
					ax.set_xticks([10,100,1000])
					ax.set_xticklabels(["$10$", "$100$", "$1000$"], fontsize=18)
					plt.xlabel("$\ell$", fontsize=20)
				if i==1:
					plt.ylabel("$\ell^2 C^{ij}_{ab}(\ell)$", fontsize=20)
				if plot_index==nbin_wl:
					plt.plot(xs, [1e-20]*len(xs), "darkviolet",linestyle="-", label="$\gamma \gamma$", lw=4.0)
					plt.plot(xg, [1e-20]*len(xg), "forestgreen",linestyle="-.", label="$\gamma \delta _g$", lw=4.0)
					plt.plot(xp, [1e-20]*len(xp), "darkred",linestyle=":", label="$\delta _g \delta _g $", lw=4.0)
					plt.legend(loc="upper left", fontsize=15)

		#plt.tight_layout()
		plt.subplots_adjust(wspace=0, hspace=0)

class fig2:
	def __init__(self, pzlist):
		self.pz_list = pzlist
		self.colours={"redmagic":"r", "skynet":"purple","skynet_biased":"forestgreen", "bpz":"slateblue", "tpz":"pink", "annz":"forestgreen", "raw_bpz":"slateblue", "uberz":"purple"}
		self.linestyles={"redmagic":"-", "skynet":"-", "skynet_biased":"--", "bpz":":", "tpz":"-", "annz":"-", "raw_bpz":":", "uberz":"-"}
		self.labels={"redmagic":"redMaGiC", "skynet":"SkyNet","skynet_biased":"SkyNet$-0.05$","annz":"ANNZ", "bpz":"BPZ", "tpz":"TPZ", "raw_bpz":"BPZ Uncalibrated", "uberz": "uberz" }



	def compare_codes(self, lsscode="redmagic", nbin=3, zmax=1.8, single_panel=False, labels=None):

		plt.subplots_adjust(wspace=0, hspace=0)

		for i in xrange(nbin+1):
			if not single_panel:
				ax=plt.subplot(nbin+1, 1, i+1)
			else:
				ax=plt.subplot(2, 1, 1)
			plt.xlim(0.0, zmax)
			print i
			if i<nbin:
				for j, pz in enumerate(self.pz_list):
					print pz.code
					if i==0:
						if not labels:
							label = self.labels[pz.code]
						else:
							label = labels[j]
					else:
						label=None
					if (pz.code!=lsscode):
						norm = np.trapz(pz.pz[i], pz.z)
						plt.plot(pz.z, pz.pz[i]/norm, self.colours[pz.code], label=label, lw=2, linestyle=self.linestyles[pz.code])
					if i==0:
						plt.legend(loc="upper right")
					plt.setp(ax.get_xticklabels(), visible=False)
					plt.setp(ax.get_yticklabels(), visible=False)
					plt.ylabel("$n(z)$", fontsize=20)


			# Do the LSS n(z) separately
			elif i==nbin:
				if single_panel:
					ax=plt.subplot(2, 1, 2)
				if not labels:
					label=self.labels[lsscode]
				else:
					label=labels[-1]
				for pz in self.pz_list:
					if pz.code==lsscode:
						z = pz.z
						p = pz.pz

				for ip, pz_to_plot in enumerate(p):
					norm = np.trapz(pz_to_plot, z)
					if ip!=0:
						label = None
					plt.plot(z, pz_to_plot/norm, self.colours[lsscode], label=label, lw=2, linestyle=self.linestyles[lsscode])
					plt.setp(ax.get_yticklabels(), visible=False)
				plt.legend(loc="upper right")


		plt.xlabel("$z$", fontsize=22)
		ax.set_xticks([0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6])
		ax.set_xticklabels(["$0.0$","$0.2$","$0.4$","$0.6$","$0.8$","$1.0$","$1.2$","$1.4$","$1.6$"], fontsize=20)
		plt.ylabel("$n(z)$", fontsize=22)


def compare_cls(datavectors, nbin_wl=3, nbin_lss=3, legend=True, additional_label="", linestyle="-"):
		import matplotlib.colors
		import matplotlib
		matplotlib.rcParams['font.family']='serif'
		matplotlib.rcParams['font.size']=16
		matplotlib.rcParams['legend.fontsize']=15
		matplotlib.rcParams['xtick.major.size'] = 10.0
		matplotlib.rcParams['ytick.major.size'] = 10.0

		plt.subplots_adjust(wspace=0, hspace=0)
		xs = datavectors[0].data["ell_shear"]
		xg = datavectors[0].data["ell_ggl"]
		xp = datavectors[0].data["ell_pos"]
		for j in xrange(1,nbin_lss+1):
			for i in xrange(1,nbin_wl+1):
				plot_index = i + (j-1)*3
				print i, j, plot_index
				ax=plt.subplot(nbin_lss, nbin_wl, plot_index )
				try:
					plt.plot(xs, datavectors[0].data["C_ee"]["bin_%d_%d"%(i,j)]/datavectors[1].data["C_ee"]["bin_%d_%d"%(i,j)], "darkviolet", linestyle=linestyle, lw=2.0)
				except:
					    print "No shear Cl"
				try:
					plt.plot(xg, datavectors[0].data["C_ne"]["bin_%d_%d"%(i,j)]/datavectors[1].data["C_ne"]["bin_%d_%d"%(i,j)], "forestgreen", linestyle=linestyle, lw=2.0)
				except:
					print "No ggl Cl"
				try:
					plt.plot(xp, datavectors[0].data["C_nn"]["bin_%d_%d"%(i,j)]/datavectors[1].data["C_nn"]["bin_%d_%d"%(i,j)], "darkred", linestyle=linestyle, lw=2.0)
				except:
					print "No position Cl"
				plt.xscale("log")
#				plt.yscale("log")
				plt.xlim(10., 2700)
#				plt.ylim(1.0e-7, 3.0)
				plt.annotate("%d,%d"%(i,j), xy=(1000.0, 6.0e-7), fontsize=18)
				#if i!=1:
					#plt.setp(ax.get_yticklabels(), visible=False)
#				elif plot_index==1:
#					ax.set_yticks([1.0e-6, 1.0e-4, 1.0e-2, 1.0])
#					ax.set_yticklabels(["$10^{-6}$", "$10^{-4}$", "$10^{-2}$", "$10^{0}$"])
#				else:
#					ax.set_yticks([1.0e-6, 1.0e-4, 1.0e-2, 1.0])
#					ax.set_yticklabels(["$10^{-6}$", "$10^{-4}$", "$10^{-2}$", "$10^{0}$"])
				if j!=nbin_lss:
					plt.setp(ax.get_xticklabels(), visible=False)
				else:
					ax.set_xticks([10,100,1000])
					ax.set_xticklabels(["$10$", "$100$", "$1000$"])
					plt.xlabel("$\ell$")
				if i==1:
					plt.ylabel("$C^{ij}_{ab, 1}(\ell)/C^{ij}_{ab, 2}(\ell)$")
				if plot_index==nbin_wl:
					if legend:
						if additional_label:
							additional_label=" "+additional_label
						#plt.plot(xs, [1.]*len(xs), "darkviolet",linestyle="-", label="Cosmic Shear"+additional_label, lw=2.0)
						#plt.plot(xg, [1.]*len(xg), "forestgreen",linestyle="-.", label="GGL"+additional_label, lw=2.0)
						#plt.plot(xp, [1.]*len(xp), "darkred",linestyle=":", label="Clust."+additional_label, lw=2.0)
						plt.legend(loc="upper left")

		#plt.tight_layout()
		plt.subplots_adjust(wspace=0.5, hspace=0.15)

class fig4:
	def extract_line(self, directname, bias=False, deltaz=True):
		print "Generating line from %s"%directname
		import glob 
		
		try:
			filename = glob.glob(directname+"/nsamp_*")[0]
			n=np.loadtxt(filename).T
		except:
			print "Could not find nsamp file"
			n=None
			filename=None
		if not bias:
			y = mc.extract_confidence_bounds(directname, deltaz=deltaz, conf="std")
		else:
			y=mc.extract_bias(directname,deltaz=deltaz)
		return n, y[0], y[1]

	def extract_bias_line(self, directname, param="s8"):
		print "Generating line from %s"%directname
		import glob 
		
		y = mc.extract_bias(directname, param=param, deltaz=True)
		return y[0], y[1]

	def composite_line(self, x, ylist, plist, ymc):
		y=[]
		for i, dzp in enumerate(x):
			# Choose two of the lines to combine
			if dzp in plist:
				sel=np.array(plist)==dzp
				y.append(ymc[sel][0])
				continue
			if dzp==x[-1]:
				y.append(ylist[-1][-1])
				continue
			elif dzp==x[0]:
				y.append(ylist[0][0])
				continue

			if dzp<plist[0]:
				y.append(ylist[0][i])
				continue
			elif dzp>plist[-1]:
				y.append(ylist[-1][i])
				continue
			upper = np.argwhere(np.array(plist)>dzp)[0][0]
			lower = np.argwhere(np.array(plist)<dzp)[-1][0]

			p1= plist[lower]
			p2 = plist[upper]

			y1 = ylist[lower][i]
			y2 = ylist[upper][i]

			print p1, p2, dzp

			d1 = (dzp-p1)
			d2 = (p2-dzp)
			dzp = abs(p2-p1)

			w2 = d1/dzp
			w1 = d2/dzp

			y.append(w1*y1+w2*y2)

		return np.array(y)


	def combine_lines(self,x,p1,y1,p2,y2):
		d1 = (np.array(x)-p1)
		d2 = (p2-np.array(x))
		dp = abs(p2-p1)

		w2 = d1/dp
		w1 = d2/dp

		y = w1*y1+w2*y2

		# If the new importance deltaz samples go outside the range between the
		# two MCMC points then revert to the importance results derived from the 
		# nearest MCMC point
		for i, val in enumerate(y): 
			if w2[i]<=0:
				y[i] = y1[i]
			if w1[i]<=0:
				y[i] = y2[i]
			else:
				continue
		return y

	def make(self, x, y, param="s8", label=True, extra_label="", probe="wl+ggl+lss", linestyle="-", lw=2.0, xlim_lower=0.0, xlim_upper=0.1, ylim_lower=0.0,ylim_upper=0.035, colour=None):
		import pylab as plt
		import matplotlib.colors
		import matplotlib
		matplotlib.rcParams['font.family']='serif'
		matplotlib.rcParams['font.size']=16

		fid=abs(fiducial[param])
		print "Normalising by fiducial value %f"%fid
		if fid==0:
			fid=1.0

		if label:
			lab=probe_labels[probe]
			if extra_label:
				lab="$%s$ "%lab+extra_label
		else:
			lab=None

		if colour:
			c=colour
		else:
			c=col[probe]

		plt.plot(x, y/fid, linestyle, color=c, label=lab, linewidth=lw, markerfacecolor="None", ms=11.0)
		plt.ylabel("Uncertainty $\sigma S_8 / S_8$ or Bias $\delta %s / %s$"%(param_names[param][1], param_names[param][1]), fontsize=20)


		plt.xlabel("Prior width $\Delta \delta z$", fontsize=20)
		plt.legend(loc="upper left", numpoints=1)
		plt.xlim(xmin=xlim_lower, xmax=xlim_upper)
		plt.ylim(ymin=ylim_lower, ymax=ylim_upper)
		plt.tight_layout()

		return 0






def generate_fig4_samples(c,number=10,xmin=0.02, xmax=0.1, write_script=True, oldprior=0.1, sgp=True, directname=""):
	"""Repeatedly apply importance sampling and resample an input MCMC chain
	   to evaluate the posterior under a range of widths of the prior on deltaz 
	"""

	deltaz=np.linspace(xmin,xmax,number)
	os.system("mkdir -p %simportance_fromdeltaz%3.2f"%(directname,oldprior))
	postprocess_command="cd /home/samuroff/cosmosis \nsource my-source \ncd - \npostprocess -o %simportance_fromdeltaz%3.2f --no-2d "%(directname,oldprior)
	N=[]
	if sgp:
		base="shear+ggl+pos"
	else:
		base="shear"
	for i in deltaz:
		c.change_deltaz_prior(i,oldprior)
		try:
			nsamp=c.write_columns(directname+base+"sig8_redsys_deltazprior%3.3f_from%3.2f_postprocessed.txt"%(i,oldprior), impweight=True)
		except:
			print "Could not process %3.4f"%i
			continue
		postprocess_command+=(directname+base+"sig8_redsys_deltazprior%3.3f_from%3.2f_postprocessed.txt "%(i,oldprior))
		N.append(nsamp)
	print "Done %f"%i

	if write_script:
		out=open(directname+"postprocess_script_fromdeltaz%3.2f"%oldprior, "w")
		out.write(postprocess_command)
		out.close()

		out=open(directname+"importance_fromdeltaz%3.2f/nsamp_fromdeltaz%3.2f"%(oldprior,oldprior), "w")
		for n in N:
			out.write("%f\n"%n)
		out.close()






class xi:
	def __init__(self, filename="/home/samuroff/cosmosis/mpp/wl_requirements/covariance/ee_y1_joachimicov.fits"):
		import fitsio as fi
		print "Getting data from :%s"%filename
		self.data={}
		infile = fitsio.FITS(filename)
		bins = np.unique(infile["xip"].read()["BIN1"])
		data = infile["xip"].read()

		self.theta = data[ (data["BIN1"]==1) & (data["BIN2"]==1) ]["ANG"]

		self.nbins = bins.max()

		for i in bins:
			for j in xrange(i, bins.max()+1):
				print "reading bin combination: %d %d"%(i,j)
				sel = (data["BIN2"]==i) & (data["BIN1"]==j)
				self.data["bin_%d_%d"%(j,i)] = data["VALUE"][sel]

		
	def make(self):
		import matplotlib.colors
		import matplotlib
		matplotlib.rcParams['font.family']='serif'
		matplotlib.rcParams['font.size']=16
		matplotlib.rcParams['legend.fontsize']=15
		matplotlib.rcParams['xtick.major.size'] = 10.0
		matplotlib.rcParams['ytick.major.size'] = 10.0

		plt.subplots_adjust(wspace=0, hspace=0)
		xs = self.theta
		for j in xrange(1,self.nbins+1):
			for i in xrange(j,self.nbins+1):
				plot_index = i + (j-1)*3
				print i, j, plot_index
				ax=plt.subplot(self.nbins, self.nbins, plot_index )
				plt.plot(xs, xs*self.data["bin_%d_%d"%(i,j)], "darkviolet", linestyle="-", lw=2.0)
				
				plt.xscale("log")
				plt.yscale("log")
				plt.xlim(0.5, 100)
				plt.ylim(1e-5,6e-4)
				plt.annotate("(%d,%d)"%(i,j), xy=(0.6, 1.8e-5), fontsize=18)
				if i!=j:
					plt.setp(ax.get_yticklabels(), visible=False)

				ax.set_yticks([1.0e-5, 1.0e-4])
				ax.set_yticklabels(["$10^{-5}$", "$10^{-4}$"], fontsize=13)
				ax.set_xticks([0.5,2,8,25,100])
				#elif plot_index==1:
					#ax.set_yticks([1.0e-6, 1.0e-4, 1.0e-2, 1.0])
					#ax.set_yticklabels(["$10^{-6}$", "$10^{-4}$", "$10^{-2}$", "$10^{0}$"], fontsize=18)
				#else:
					#ax.set_yticks([2.0e-6, 1.0e-4, 1.0e-2, 1.0])
					#ax.set_yticklabels(["$10^{-6}$", "$10^{-4}$", "$10^{-2}$", "$10^{0}$"], fontsize=18)
				if j!=i:
					plt.setp(ax.get_xticklabels(), visible=False)
				else:
					if i==j==self.nbins:
						ax.set_xticklabels(["$0.5$", "$2.0$", "$8.0$","$25.0$", "$100.0$"], fontsize=13)
						plt.xlabel(r"angular scale $\theta$ / arcmin", fontsize=14)
					else:
						ax.set_xticklabels(["$0.5$", "$2.0$", "$8.0$","$25.0$"], fontsize=13)
					
				if i==j:
					plt.ylabel(r"$\theta \xi_+ ^{i,j} (\theta)$", fontsize=14)

		#plt.tight_layout()
		plt.subplots_adjust(wspace=0, hspace=0)











labels={'snr' : ("$SNR$", "snr"), "rgpp":("$R_{gpp}/R_{p}$", "mean_rgpp_rp"), "e1": ("$e_1$", "e1"), "e2": ("$e_2$", "e2"), "dphi": (r"shape-shape misalignment $\Delta \phi$ / $\pi$ rad", "dphi"), "dbeta": (r"shape-direction misalignment $\Delta \beta$ / $\pi$ rad", "dbeta"), "m":("$m \equiv (m_1+m_2)/2$", "m"), "m11":("$m_{1}$", "m11"), "m22":("$m_{2}$", "m22"), "c":("$c \equiv (c_1+c_2)/2$", "c"), "c11":("$c_{1}$", "c11"), "c22":("$c_{2}$", "c22"), "alpha":(r"$\alpha$", "alpha")  }

class im3shape_results_plots:
	def obs_vs_obs(self, name1, name2, nbins, binning="equal_number", xlabel=None, ylabel=None, label=None, ls="-", xlim=None, ylim=None, fmt="o", colour="purple", scatter=False, mean=True, return_vals=False,logscale=False,refline=None, ydata=None, xdata=None, plot=True):
		import pylab as plt
		col1=None
		lab1=None
		col2=None
		lab2=None
		xmin=None
		xmax=None
		ymin=None
		ymax=None

		# Do the labelling
		try:
			lab1,col1=labels[name1]
		except:
			print "No label for %s found"%name1
			lab1,col1=name1,name1
		try:
			lab2,col2=labels[name2]
		except:
			print "No label for %s found"%name2
			lab2,col2=name2,name2

		if xlabel is not None: lab1 = xlabel
		if ylabel is not None: lab2 = ylabel

		if plot:
			plt.xlabel(lab1)
			plt.ylabel(lab2)

		# Deal with the axis bounds
		if xlim is not None:
			xmin=xlim[0]
			xmax=xlim[1]
			if plot:
				plt.xlim(xmin=xmin,xmax=xmax)
		if ylim is not None:
			ymin=ylim[0]
			ymax=ylim[1]
			if plot:
				plt.ylim(ymin=ymin,ymax=ymax)

		
		data=self.res
		if col1 not in data.dtype.names:
			data = arr.add_col(data, col1, self.truth[col1])
		if ydata is None:
			if col2 not in data.dtype.names:
				data = arr.add_col(data, col2, self.truth[col2])
		else:
			data= arr.add_col(data, col2, ydata)

		if xmin is not None:
			sel1 = (self.res[col1]>xmin)
		else:
			sel1 = np.ones_like(self.res).astype(bool)
		if xmax is not None:
			sel2 = (self.res[col1]<xmax)
		else:
			sel2 = np.ones_like(self.res).astype(bool)
		if ymin is not None:
			sel3 = (self.res[col2]>ymin)
		else:
			sel3 = np.ones_like(self.res).astype(bool)
		if ymax is not None:
			sel4 = (self.res[col2]<ymax)
		else:
			sel4 = np.ones_like(self.res).astype(bool)

		try:
			sel5 = np.isfinite(self.res[col2])
		except:
			sel5 = np.isfinite(self.res["e1"])

		data = data[sel1 & sel2 & sel3 & sel4 & sel5]
		print "Total galaxies:%d"%data.size

		if scatter:
			plt.scatter(data[col1], data[col2], color=colour)

		elif mean:
			if binning=="equal_number":
				bin_edges=di.find_bin_edges(data[col1], nbins)
			else:
				bin_edges=np.array(binning)
			x=(bin_edges[1:]+bin_edges[:-1])/2.
			y=[]
			err=[]
			n=[]

			for i,bounds in enumerate(zip(bin_edges[:-1],bin_edges[1:])):
				upper = bounds[1]
				lower = bounds[0]
				print "bin %d [%f,%f]"%((i+1),lower,upper)
				sel = (data[col1]>lower) & (data[col1]<upper)
				if (name2 is not "m1") and (name2 is not "m2") and (name2 is not "m") and (name2 is not "c1") and (name2 is not "c2") and (name2 is not "alpha") and (name2 is not "alpha1") and (name2 is not "alpha2"):
					q = data[col2][sel]
					print q.size
					n.append(q.size)
					e = np.std(q)/np.sqrt(q.size)
				else:
					tr = self.truth[sel1 & sel2 & sel3 & sel4 & sel5]
					try: 
						q,e= di.truth_plots(data[sel], tr[sel], nbins=5, mode="none", psf_leakage=True, bin='', return_vals=name2, use_intrinsic=False, equal_number_bins=True)
					except:
						q,e=0.0,0.0
					print q
				
				y.append(np.mean(q))
				err.append(e)

			if plot:
				plt.errorbar(x,y,yerr=err,lw=2.5,color=colour,fmt=fmt)
				plt.plot(x,y,lw=2.5,color=colour,ls=ls,label=label)
				if logscale:
					plt.xscale("log")
				if refline is not None:
					plt.axhline(refline,color="k")

		if return_vals:
			return np.array(x),np.array(y),np.array(err),np.array(n)
		else:
			return 0

	def truth_vs_obs(self, name1, name2, nbins, binning="equal_number", label=None, ls="-", xlim=None, ylim=None, fmt="o", colour="purple", scatter=False, mean=True, return_vals=False,logscale=False,refline=None):
		import pylab as plt
		col1=None
		lab1=None
		col2=None
		lab2=None
		xmin=None
		xmax=None
		ymin=None
		ymax=None

		# Do the labelling
		try:
			lab1,col1=labels[name1]
		except:
			print "No label for %s found"%name1
			lab1,col1=name1,name1
		try:
			lab2,col2=labels[name2]
		except:
			print "No label for %s found"%name2
			lab2,col2=name2,name2

		plt.xlabel(lab1)
		plt.ylabel(lab2)

		# Deal with the axis bounds
		if xlim is not None:
			xmin=xlim[0]
			xmax=xlim[1]
			plt.xlim(xmin=xmin,xmax=xmax)
		if ylim is not None:
			ymin=ylim[0]
			ymax=ylim[1]
			plt.ylim(ymin=ymin,ymax=ymax)

		if xmin is not None:
			sel1 = (self.truth[col1]>xmin)
		else:
			sel1 = np.ones_like(self.truth).astype(bool)
		if xmax is not None:
			sel2 = (self.truth[col1]<xmax)
		else:
			sel2 = np.ones_like(self.res).astype(bool)
		if ymin is not None:
			sel3 = (self.res[col2]>ymin)
		else:
			sel3 = np.ones_like(self.res).astype(bool)
		if ymax is not None:
			sel4 = (self.res[col2]<ymax)
		else:
			sel4 = np.ones_like(self.res).astype(bool)

		try:
			sel5 = np.isfinite(self.res[col2])
		except:
			sel5 = np.isfinite(self.res["e1"])

		datax = self.truth[sel1 & sel2 & sel3 & sel4 & sel5]
		datay = self.res[sel1 & sel2 & sel3 & sel4 & sel5]
		print "Total galaxies:%d"%datay.size

		if binning=="equal_number":
			bin_edges=di.find_bin_edges(datax[col1], nbins)
		else:
			bin_edges = np.array(binning)

		if bin_edges.size!=nbins:
			raise ValueError("You asked for %d bins but gave %d bin edges"%(nbins, bin_edges.size))

		x=(bin_edges[1:]+bin_edges[:-1])/2.
		y=[]
		err=[]
		n=[]

		for i,bounds in enumerate(zip(bin_edges[:-1],bin_edges[1:])):
			upper = bounds[1]
			lower = bounds[0]
			print "bin %d [%f,%f]"%((i+1),lower,upper)
			sel = (datax[col1]>lower) & (datax[col1]<upper)
			if (name2 is not "m1") and (name2 is not "m2") and (name2 is not "m"):
				q = datay[col2][sel]
				print q.size
				n.append(q.size)
				e = np.std(q)/np.sqrt(q.size)
			else:
				tr = self.truth[sel1 & sel2 & sel3 & sel4 & sel5]
				try:
					(m11,c11,em11,ec11),(m22,c22,em22,ec22),(m12,c12,em12,ec12),(m21,c21,em21,ec21),(a11,c11psf,ea11,ecpsf11),(a22,cpsf22,ea22,ecpsf22),fb = di.truth_plots(datay[sel], tr[sel], nbins=5, mode="none", psf_leakage=True, bin='', return_vals=True, use_intrinsic=False, equal_number_bins=True)
				except:
					(m11,c11,em11,ec11),(m22,c22,em22,ec22),(m12,c12,em12,ec12),(m21,c21,em21,ec21),(a11,c11psf,ea11,ecpsf11),(a22,cpsf22,ea22,ecpsf22),fb = (0.0,0.0,0.0,0.0),(0.0,0.0,0.0,0.0),(0.0,0.0,0.0,0.0),(0.0,0.0,0.0,0.0),(0.0,0.0,0.0,0.0),(0.0,0.0,0.0,0.0),0.0
				if name2=="m":
					q = np.array([(m11+m22)/2.0])
					e = np.sqrt(em11*em11+em22*em22)
				elif name2=="m1":
					q = np.array([m11])
					e = em11
				elif name2=="m2":
					q = np.array([m22])
					e = em22
				print q
			
			y.append(np.mean(q))
			err.append(e)

		plt.errorbar(x,y,yerr=err,lw=2.5,color=colour,fmt=fmt)
		plt.plot(x,y,lw=2.5,color=colour, label=label, ls=ls)
		if logscale:
			plt.xscale("log")
		if refline is not None:
			plt.axhline(refline,color="k")

		if return_vals:
			return np.array(x),np.array(y),np.array(err),np.array(n)
		else:
			return 0

	def bias_vs_obs(self, name, bias, nbins, binning="equal_number", ellipticity_name="e", error_type="bootstrap", xlabel=None, label=None, ls="-", xlim=None, ylim=None, fmt="o", colour="purple", return_vals=False, logscale=False, refline=0, extra_selection=[], xdata=None):
		import pylab as plt
		col1=None
		lab1=None
		xmin=None
		xmax=None
		ymin=None
		ymax=None

		# Do the labelling
		lab2,col2 = labels[bias]
		try:
			lab1,col1=labels[name]
		except:
			print "No label for %s found"%name
			lab1,col1=name,name

		if xlabel is not None:
			lab1 = xlabel

		plt.xlabel(lab1)
		plt.ylabel(lab2)

		# Deal with the axis bounds
		axis_lower, axis_upper = axis_bounds(xlim,ylim,self.res[name])

		if len(extra_selection)>0:
			sel = extra_selection
		else:
			sel = np.ones_like(self.res["e1"]).astype(bool)

		data = self.res[axis_lower & axis_upper & sel]
		if hasattr(self, "truth"):
			tr = self.truth[axis_lower & axis_upper & sel]
			data = arr.add_col(data,"true_g1",tr["true_g1"])
			data = arr.add_col(data,"true_g2",tr["true_g2"])
			data = arr.add_col(data,"intrinsic_sheared_e1",tr["intrinsic_e1"]+tr["true_g1"])
			data = arr.add_col(data,"intrinsic_sheared_e2",tr["intrinsic_e2"]+tr["true_g2"])

		if xdata is not None:
			data[name] = xdata


		print "Total galaxies:%d"%data.size

		#Now the binning
		if binning=="equal_number":
			bin_edges=di.find_bin_edges(data[col1], nbins)
		else:
			bin_edges=np.array(binning)
		x=(bin_edges[1:]+bin_edges[:-1])/2.
		y=[]
		err=[]
		n=[]

		for i,bounds in enumerate(zip(bin_edges[:-1],bin_edges[1:])):
			upper = bounds[1]
			lower = bounds[0]
			print "bin %d [%f,%f]"%((i+1),lower,upper)
			bin = (data[col1]>lower) & (data[col1]<upper)
			# Case where we are calculating biases
			if bias not in data.dtype.names:
				# Do the linear fits
				# Should return a dictionary, each element of which is a value then an uncertainty
				b = di.get_bias(data[bin], nbins=5, ellipticity_name=ellipticity_name, binning="equal_number", names=["m","c","m11","m22","c11","c22"])
				# Repeat them if the errorbars need to come from bootstrapping
				if error_type=="bootstrap":
					try:
						error = di.bootstrap_error(6, data[bin], di.get_bias, additional_args=["names", "silent", "ellipticity_name"], additional_argvals=[bias, True, ellipticity_name])
					except:
						error = -1
				else:
					error = b[bias][1]

			# Case where we have precomputed biases
			else:
				b = {bias:( np.mean(data[bias]), np.std(data[bias])/ (data[bias].size**0.5) ) }
					
			y.append(b[bias][0])
			err.append(error)
			n.append(data[bin].size)

		print x
		print y
		print err

		plt.errorbar(x,y,yerr=err,lw=2.5,color=colour,fmt=fmt)
		plt.plot(x,y,lw=2.5,color=colour,ls=ls,label=label)
		if logscale:
			plt.xscale("log")
		if refline is not None:
			plt.axhline(refline,color="k")

		if return_vals:
			return np.array(x),np.array(y),np.array(err),np.array(n)
		else:
			return 0

	def alpha_vs_obs(self, name, bias, nbins, binning="equal_number", ellipticity_name="e", error_type="bootstrap", xlabel=None, label=None, ls="-", xlim=None, ylim=None, fmt="o", colour="purple", return_vals=False, logscale=False, refline=0):
		import pylab as plt
		col1=None
		lab1=None
		xmin=None
		xmax=None
		ymin=None
		ymax=None

		# Do the labelling
		try:
			lab1,col1=labels[name]
		except:
			print "No label for %s found"%name
			lab1,col1=name,name
		lab2,col2=labels["m"]

		if xlabel is not None:
			lab1 = xlabel

		plt.xlabel(lab1)
		plt.ylabel(lab2)

		# Deal with the axis bounds
		axis_lower, axis_upper = axis_bounds(xlim,ylim,self.res[name])

		data = self.res[axis_lower & axis_upper]
		if hasattr(self, "truth"):
			tr = self.truth[axis_lower & axis_upper]
			data = arr.add_col(data,"true_g1",tr["true_g1"])
			data = arr.add_col(data,"true_g2",tr["true_g2"])
			data = arr.add_col(data,"intrinsic_sheared_e1",tr["intrinsic_e1"]+tr["true_g1"])
			data = arr.add_col(data,"intrinsic_sheared_e2",tr["intrinsic_e2"]+tr["true_g2"])


		print "Total galaxies:%d"%data.size

		#Now the binning
		if binning=="equal_number":
			bin_edges=di.find_bin_edges(data[col1], nbins)
		else:
			bin_edges=np.array(binning)
		x=(bin_edges[1:]+bin_edges[:-1])/2.
		y=[]
		err=[]
		n=[]

		for i,bounds in enumerate(zip(bin_edges[:-1],bin_edges[1:])):
			upper = bounds[1]
			lower = bounds[0]
			print "bin %d [%f,%f]"%((i+1),lower,upper)
			bin = (data[col1]>lower) & (data[col1]<upper)
			# Case where we are calculating biases
			if bias not in data.dtype.names:
				# Do the linear fits
				# Should return a dictionary, each element of which is a value then an uncertainty
				b = di.get_alpha(data[bin], nbins=5, ellipticity_name=ellipticity_name, binning="equal_number", names=["alpha","c","alpha11","alpha22","c11","c22"], use_weights=False)
				# Repeat them if the errorbars need to come from bootstrapping
				if error_type=="bootstrap":
					error = di.bootstrap_error(6, data[bin], di.get_alpha, additional_args=["names", "silent", "ellipticity_name"], additional_argvals=[bias, True, ellipticity_name])
				else:
					error = b[bias][1]

			# Case where we have precomputed biases
			else:
				b = {bias:( np.mean(data[bias]), np.std(data[bias])/ (data[bias].size**0.5) ) }
					
			y.append(b[bias][0])
			err.append(error)

		print x
		print y
		print err

		plt.errorbar(x,y,yerr=err,lw=2.5,color=colour,fmt=fmt)
		plt.plot(x,y,lw=2.5,color=colour,ls=ls,label=label)
		if logscale:
			plt.xscale("log")
		if refline is not None:
			plt.axhline(refline,color="k")

		if return_vals:
			return np.array(x),np.array(y),np.array(err),np.array(n)
		else:
			return 0

	def errorbar_vs_ngal(self, bias, nbins, ellipticity_name="e", error_type="bootstrap", xlabel=None, label=None, ls="-", colour="purple", return_vals=False, logscale=False, refline=0, show=True):
		import pylab as plt
		col1=None
		lab1=None
		xmin=None
		xmax=None
		ymin=None
		ymax=None

		bias_base=bias[0]
		bias_function_lookup={"a": di.get_alpha, "m": di.get_bias, "c": di.get_bias}
		bias_function = bias_function_lookup[bias_base]

		# Do the labelling
		lab2,col2=labels[bias]

		data = self.res
		if hasattr(self, "truth"):
			tr = self.truth
			data = arr.add_col(data,"true_g1",tr["true_g1"])
			data = arr.add_col(data,"true_g2",tr["true_g2"])
			data = arr.add_col(data,"intrinsic_sheared_e1",tr["intrinsic_e1"]+tr["true_g1"])
			data = arr.add_col(data,"intrinsic_sheared_e2",tr["intrinsic_e2"]+tr["true_g2"])


		print "Total galaxies:%d"%data.size

		#Now the binning
		x = np.linspace(0,1,nbins+1)[1:]
		y=[]
		err=[]

		for i, frac in enumerate(x):

			ngal = frac*data.size
			subset = np.random.choice(data,ngal.astype(int))
			# Case where we are calculating biases
			if bias not in data.dtype.names:
				# Do the linear fits
				# Should return a dictionary, each element of which is a value then an uncertainty
				b = bias_function(subset, nbins=5, ellipticity_name=ellipticity_name, silent=np.invert(show), binning="equal_number", names=[bias])
				# Repeat them if the errorbars need to come from bootstrapping
				if error_type=="bootstrap":
					error = di.bootstrap_error(6, subset, bias_function, additional_args=["names", "silent", "ellipticity_name"], additional_argvals=[bias, True, ellipticity_name])
				else:
					error = b[bias][1]

			# Case where we have precomputed biases
			else:
				b = {bias:( np.mean(data[bias]), np.std(data[bias])/ (data[bias].size**0.5) ) }
					
			y.append(error)
			#err.append(error)

		x = (np.array(x).astype(float) * data.size)/1e6
		y = np.array(y)

		if show:
			plt.plot(x,y,lw=2.5,color=colour,ls=ls,label=label)
			plt.xlim(0,x.max())
			plt.ylim(0,y.max())
			if logscale:
				plt.xscale("log")
			if refline is not None:
				plt.axhline(refline,color="k")

			plt.ylabel(r"statistical error $\sigma \lbrack %s \rbrack$"%lab2.replace("$",""))
			plt.xlabel("number of objects $N_{gal}$ (M galaxies)")

		if return_vals:
			return np.array(x),np.array(y)
		else:
			return 0

	def dot_plot(self, nbins, data=None, weights=True):
		if data is None:
			data = shapecat(res="/share/des/disc2/y1/y1a1-im3shape-r-1-1-1.fits")
			data.load(truth=False)
			data.apply_infocuts()

		sbins=np.logspace(data.res["snr"].min(),data.res["snr"].max() , nbins+1)
		rbins=np.linspace(data.res["mean_rgpp_rp"].min(),data.res["mean_rgpp_rp"].max(), nbins+1)

		# Grid the weights
		# This is potentially slow- use weights=False to skip this
		if weights:
			wts=np.zeros((nbins,nbins))
			for i in xrange(nbins):
				for j in xrange(nbins):
					print i,j
					slower=sbins[i] ; supper=sbins[i+1]
					rlower=rbins[j] ; rupper=rbins[j+1]
					sel1 = (data.res["snr"]>slower) & (data.res["snr"]<supper)
					sel2 = (data.res["mean_rgpp_rp"]>rlower) & (data.res["mean_rgpp_rp"]<rupper)

					w = dat.res["weight"][sel1 & sel2]

					if w.size==0:
						wts[i,j]=0
					else:
						wts[i,j] = np.mean()

		else:
			wts = np.ones((nbins,nbins))

		x = np.sqrt(sbins[1:]*sbins[:-1])
		y = (rbins[1:]+rbins[:-1])/2.0

		simgrid = np.histogram2d(sim.res["snr"], sim.res["mean_rgpp_rp"], bins=[sbins,rbins])
		datgrid = np.histogram2d(dat.res["snr"], dat.res["mean_rgpp_rp"], bins=[sbins,rbins])
		k = (1.0*data.res.size)/self.res.size

		xy = np.meshgrid(x,y)

		snr = xy[0].flatten()
		rgpp = xy[1].flatten()

		plt.scatter(np.log10(snr),rgpp, c=k*simgrid[0].flatten()/datgrid[0].flatten(), s=(20*wts.flatten()))

	








	def whisker_plot(self, nobj=None, thin=None, etype="galaxies", positions="world", colour="purple", lw=0.5, label=None, newplot=1, zoom=None, unit_length=None):

		plt.figure(newplot)

		if etype=="galaxies":
			z1 = self.res["e1"]
			z2 = self.res["e2"]
		elif etype=="psf":
			z1 = self.res["mean_psf_e1_sky"]
			z2 = self.res["mean_psf_e2_sky"]
		elif etype=="shear":
			z1 = self.truth["true_g1"]
			z2 = self.truth["true_g2"]

		if nobj is None:
			nobj=self.res.size
			
		x0,y0 = self.get_positions(positions)
		x,y = x0,y0

		if nobj is not None:
			print "using first %d objects"%nobj
			x=x[:nobj]
			y=y[:nobj]
		if thin is not None:
			x=x[::thin]
			y=y[::thin]
			nobj = x.size
			print "applying thinning factor of %d."%thin

		print "Will make whisker plot using %d objects"%nobj

		if unit_length is None:
			if etype=="galaxies":
				unit_length_x = 180.0
				unit_length_y = 180.0
			elif etype=="psf":
				unit_length_x = 4800.0
				unit_length_y = 4800.0
			elif etype=="shear":
				unit_length_x = 4800.0
				unit_length_y = 4800.0
		else:
			unit_length_x = unit_length
			unit_length_y = unit_length

		for i, cent in enumerate(zip(x,y)):
			e1 = z1[i]
			e2 = z2[i]

			if e1==-9999:
				e1 = 0
			if e2==-9999:
				e2 = 0

			xcent,ycent = cent

			a = [xcent-(e1*unit_length_x), xcent+(e1*unit_length_x)]
			b = [ycent-(e2*unit_length_y), ycent+(e2*unit_length_y)]

			plt.plot(a,b,lw=lw, color=colour)

		if zoom is None:
			zoom=1.0
		if not isinstance(zoom,float):
			plt.xlim(xmin=zoom[0],xmax=zoom[1])
			plt.ylim(ymin=zoom[0],ymax=zoom[1])
		else:
			plt.xlim(xmin=x0.min(),xmax=x0.max()/zoom)
			plt.ylim(ymin=y0.min(),ymax=y0.max()/zoom)

		plt.title(label)

		return 0

	def neighbour_cartoon(self, ncat, iobj, angles=True):

		fig = plt.figure()
		ax = fig.add_subplot(111, aspect='equal')

		cent = self.res[iobj]
		ul=18
		f= 60*60/0.27
		boxsize = self.res[iobj]["stamp_size"]

		x_c = self.truth["ra"][iobj]*f
		y_c = self.truth["dec"][iobj]*f
		e1_c = self.res["e1"][iobj]
		e2_c = self.res["e2"][iobj]
		phi_c = np.angle(e1_c+1j*e2_c)/2
		x_n = ncat.truth["ra"][iobj]*f
		y_n = ncat.truth["dec"][iobj]*f
		e1_n = ncat.res["e1"][iobj]
		e2_n = ncat.res["e2"][iobj]
		phi_n = np.angle(e1_n+1j*e2_n)/2

		e_c = di.quad(e1_c, e2_c)
		e_n = di.quad(e1_n, e2_n)

		ex_c= e_c*np.cos(phi_c) * ul
		ey_c= e_c*np.sin(phi_c) * ul
		ex_n= e_n*np.cos(phi_n) * ul
		ey_n= e_n*np.sin(phi_n) * ul

		plt.plot([x_c - ex_c, x_c + ex_c ], [y_c-ey_c, y_c+ey_c ], "-", color="steelblue", lw=2.0)
		plt.plot([x_n - ex_n, x_n + ex_n ], [y_n-ey_n, y_n+ey_n ], "-", color="forestgreen", lw=2.0)
		plt.plot([x_c - ex_n, x_c + ex_n ], [y_c-ey_n, y_c+ey_n ], ":", color="forestgreen", lw=2.0)
		plt.plot([x_c, x_n], [y_c, y_n], "-", color="purple", lw=2.5)
		from matplotlib.patches import Ellipse

		gal=Ellipse(xy=[x_c,y_c], width=2*ex_c, height=2*ey_c, angle=0)
		plt.rtist(gal)
		gal.set_facecolor("none")


		#Mark on the stamp edges
		ax.plot([x_c-boxsize/2, x_c+boxsize/2], [y_c+boxsize/2, y_c+boxsize/2], "k--") # Upper
		ax.plot([x_c-boxsize/2, x_c+boxsize/2], [y_c-boxsize/2, y_c-boxsize/2], "k--") # Lower
		ax.plot([x_c-boxsize/2, x_c-boxsize/2], [y_c-boxsize/2, y_c+boxsize/2], "k--") # Left
		plt.plot([x_c+boxsize/2, x_c+boxsize/2], [y_c-boxsize/2, y_c+boxsize/2], "k--") # Right

		neigh=Ellipse(xy=[x_n,y_n], width=2*ex_n, height=2*ey_n, angle=0)
		plt.add_artist(neigh)
		neigh.set_facecolor("none")

		if angles:
			beta=self.res["dbeta"][iobj]
			phi=self.res["dphi"][iobj]
			plt.title(r"$\Delta \phi = %1.3f \pi$ rad, $\Delta \beta = %1.3f \pi$ rad"%(phi,beta))

		plt.gca().set_aspect('equal', adjustable='box')
		plt.draw()

	def bias_fit_vs_pts(self, table=None, bias_name="m", output=None, use_rbf=False, do_half=0):
		import tools.nbc as cal
		cal.show_table(table, legend=True, name=bias_name*(bias_name=="m") + "alpha"*(bias_name=="a"), do_half=do_half)

		bt = fitsio.FITS(table)[1].read()
		snr = np.linspace(np.unique(bt["snr_lower"])[0], 250, 1000)
		rgpp = (np.unique(bt["rgp_lower"])+np.unique(bt["rgp_upper"]))/2
		nbins = rgpp.size


		plt.xscale("log")
		colours=["purple", "forestgreen", "steelblue", "pink", "darkred", "midnightblue", "gray", "sienna", "olive", "darkviolet", "cyan"]
		for i,r in enumerate(rgpp):
			if do_half==1 and i>nbins/2:
				continue
			elif do_half==2 and i<nbins/2:
				continue
			x = np.array([snr,np.array([r]*snr.size) ]).T
			if not use_rbf:
				com2 = "a0 "+", a%d "*self.optimised_coefficients_m[1:].size +"=tuple(self.optimised_coefficients_%s)"%bias_name
				com2 = com2%tuple(np.linspace(1,17,17))
				exec com2

				com = "cal.eval_%s(x"%bias_name + ", a%d "*self.optimised_coefficients_m.size+")"
				com = com%tuple(np.linspace(0,17,18))
				exec "bias=%s"%com
			else:
				bias = self.do_rbf_interpolation(bias_name,snr,np.array([r]*snr.size))

			plt.plot(snr+np.log(i+1), bias, color=colours[i])

		if bias_name is "m":
			plt.ylim(-0.5,0.18)
			plt.legend(loc="lower right")
		else:
			plt.ylim(-0.5,1.5)
			plt.legend(loc="upper right")

		plt.savefig(output)
		plt.close()

	def redshift_diagnostic(self, bias="m", split_half=2, colour="purple", fmt=["o","D"], ls="none", label=None, ellipticity_name="e", apply_calibration=False, error_type="bootstrap", nbins=5, legend=True):
		bins = np.linspace(0,1.35,nbins+1)
		lower = bins[:-1]
		upper = bins[1:]

		if "m" in bias:
			xmin,xmax= -1,1
		else:
			xmin,xmax = -0.015, 0.020, 
		exec "data = self.res%d"%split_half
		exec "truth = self.truth%d"%split_half

		vec1 = []
		err_vec1 = []
		vec2 = []
		err_vec2 = []
		z = []

		for i, edges in enumerate(zip(lower,upper)):
			sel = (truth["cosmos_photoz"]>edges[0]) & (truth["cosmos_photoz"]<edges[1])
			if bias in ["alpha", "alpha11", "alpha22"]:
				bias_function = di.get_alpha
				xdata= data
				names = ["alpha","c","alpha11","alpha22","c11","c22"]
			else:
				bias_function = di.get_bias
				xdata = truth
				names = ["m","c","m11","m22","c11","c22"]
			
			b = bias_function(xdata[sel], data[sel], nbins=5, apply_calibration=apply_calibration, ellipticity_name=ellipticity_name, binning="equal_number", names=names)
			# Repeat them if the errorbars need to come from bootstrapping
			if error_type=="bootstrap":
				error1 = di.bootstrap_error(6, (xdata[sel], data[sel]), bias_function, additional_args=["names", "nbins", "apply_calibration", "silent", "ellipticity_name", "xlim"], additional_argvals=[names[2], 6, apply_calibration, True, ellipticity_name, (xmin,xmax)])
				error2 = di.bootstrap_error(6, (xdata[sel], data[sel]), bias_function, additional_args=["names", "nbins", "apply_calibration", "silent", "ellipticity_name", "xlim"], additional_argvals=[names[3], 6, apply_calibration, True, ellipticity_name, (xmin,xmax)])
			else:
				error1 = b["%s11"%bias][1]
				error2 = b["%s22"%bias][1]

			z.append(truth[sel]["cosmos_photoz"].mean())
			vec1.append(b["%s11"%bias][0])
			err_vec1.append(error1)
			vec2.append(b["%s22"%bias][0])
			err_vec2.append(error2)

		plt.errorbar(z,vec1,err_vec1, lw=2.5, ls=ls, color=colour, label="$%s_1$ %s"%(r'\alpha'*(bias=="alpha")+"m"*(bias!="alpha"),label), fmt=fmt[0] )
		plt.errorbar(z,vec2,err_vec2, lw=2.5, ls=ls, color=colour, label="$%s_2$ %s"%(r'\alpha'*(bias=="alpha")+"m"*(bias!="alpha"),label), fmt=fmt[1] )
		if bias=="m":
			plt.axhspan(-0.02,0.02,color="forestgreen", alpha=0.3)
		plt.axhline(0, color="k", lw=2.5)
		#plt.ylim(-0.7,0.1)
		plt.xlabel("Redshift $z$")
		return z,vec1,err_vec1, vec2,err_vec2


def axis_bounds(xlim,ylim,xdata):
	if xlim is not None:
		xmin=xlim[0]
		xmax=xlim[1]
		plt.xlim(xmin=xmin,xmax=xmax)
	else:
		xmin=None
		xmax=None
	
	if ylim is not None:
		ymin=ylim[0]
		ymax=ylim[1]
		plt.ylim(ymin=ymin,ymax=ymax)
	else:
		ymin=None
		ymax=None

	if xmin is not None:
		sel1 = (xdata>xmin)
	else:
		sel1 = np.ones_like(xdata).astype(bool)
	if xmax is not None:
		sel2 = (xdata<xmax)
	else:
		sel2 = np.ones_like(xdata).astype(bool)

	return sel1, sel2



import matplotlib as mpl
def reverse_colourmap(cmap, name = 'my_cmap_r'):

    reverse = []
    k = []   

    for key in cmap._segmentdata:    
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:                    
            data.append((1-t[0],t[2],t[1]))            
        reverse.append(sorted(data))    

    LinearL = dict(zip(k,reverse))
    my_cmap_r = mpl.colors.LinearSegmentedColormap(name, LinearL) 
    return my_cmap_r

def histograms(names, data, outdir="/home/samuroff/shear_pipeline/end-to-end/plots/v2.2/new", data2=None, thin=None):

	os.system("mkdir -p %s"%outdir)

	print "Making plots:"

	if "snr" in names:
		print "-- SNR"
		plt.hist(plt.log10(data["snr"]), histtype="step", bins=np.linspace(0.8,3.0, 50), normed=1, lw=2.5, color="purple", label="v2 Sim")
		if data2 is not None:
			plt.hist(plt.log10(data2["snr"]), histtype="step", bins=np.linspace(0.8,3.0,50), normed=1, lw=2.5, color="steelblue", label="v02 Data")
			plt.legend(loc="upper right")
		plt.xlabel("Signal-to-Noise $\log (SNR_w)$")
		plt.savefig("%s/snr-hist-v2sim-vs-y1v2data.png"%outdir)
		plt.close()

	if "rgpp" in names:
		print "-- Rgpp"
		plt.hist(data["mean_rgpp_rp"], histtype="step", bins=np.linspace(0.6,2.5, 50), normed=1, lw=2.5, color="purple", label="v2 Sim")
		if data2 is not None:
			plt.hist(data2["mean_rgpp_rp"], histtype="step", bins=np.linspace(0.6,2.5,50), normed=1, lw=2.5, color="steelblue", label="v02 Data")
			plt.legend(loc="upper right")
		plt.xlabel("Size $R_{gpp} / R_p $")
		plt.savefig("%s/rgpp_rp-hist-v2sim-vs-y1v2data.png"%outdir)
		plt.close()

	if "flux" in names:
		print "-- Mean Flux"
		plt.hist(data["mean_flux"], histtype="step", bins=np.linspace(0.,3000, 50), normed=1, lw=2.5, color="purple", label="v2 Sim")
		if data2 is not None:
			plt.hist(data2["mean_flux"], histtype="step", bins=np.linspace(0.,3000,50), normed=1, lw=2.5, color="steelblue", label="v02 Data")
			plt.legend(loc="upper right")
		plt.xlabel("Mean Flux $f$")
		plt.savefig("%s/flux-hist-v2sim-vs-y1v2data.png"%outdir)
		plt.close()

	if "e" in names:
		print "-- Ellipticity"
		fig = plt.figure(1)
		ax = fig.add_subplot(131)
		ax.set_xticks([-0.5,0,0.5,1.0])
		plt.setp(ax.get_yticklabels(), visible=False)
		ax.hist(data["e1"], histtype="step", bins=np.linspace(-1,1, 50), normed=1, lw=2.5, color="purple", label="v2 Sim")
		if data2 is not None:
			ax.hist(data2["e1"], histtype="step", bins=np.linspace(-1,1,50), normed=1, lw=2.5, color="steelblue", label="v02 Data")
		ax.set_xlabel("Ellipticity $e_1$")

		ax2 = fig.add_subplot(132)
		plt.setp(ax2.get_yticklabels(), visible=False)
		ax2.set_xticks([-0.5,0,0.5,1.0])
		ax2.hist(data["e2"], histtype="step", bins=np.linspace(-1,1, 50), normed=1, lw=2.5, color="purple", label="v2 Sim")
		if data2 is not None:
			ax2.hist(data2["e2"], histtype="step", bins=np.linspace(-1,1,50), normed=1, lw=2.5, color="steelblue", label="v02 Data")
		ax2.set_xlabel("Ellipticity $e_2$")

		ax3 = fig.add_subplot(133)
		plt.setp(ax3.get_yticklabels(), visible=False)
		ax3.set_xticks([0.25,0.5,0.75,1.0])
		e = np.sqrt(data["e1"]*data["e1"] + data["e2"]*data["e2"])
		ax3.hist(e/(e.max()), histtype="step", bins=np.linspace(0.,1, 50), normed=1, lw=2.5, color="purple", label="v2 Sim")
		if data2 is not None:
			e = np.sqrt(data2["e1"]*data2["e1"] + data2["e2"]*data2["e2"])
			ax3.hist(e/e.max(), histtype="step", bins=np.linspace(0.,1,50), normed=1, lw=2.5, color="steelblue", label="v02 Data")
			matplotlib.rcParams['legend.fontsize']=10
			ax3.legend(loc="upper right")
		ax3.set_xlabel("Ellipticity $|e|$")

		plt.subplots_adjust(hspace=0, wspace=0, left=0.01,right=0.97, top=0.62)
		plt.savefig("%s/ellipticity-hist-v2sim-vs-y1v2data.png"%outdir)
		plt.close()

	if "psfe" in names:
		print "-- PSF Ellipticity"
		plt.subplots_adjust(wspace=0, hspace=0)
		fig = plt.figure(2)
		ax = fig.add_subplot(131)
		plt.setp(ax.get_yticklabels(), visible=False)
		ax.set_xticks([-0.05,0,0.025,0.05])
		ax.hist(data["mean_psf_e1_sky"], histtype="step", bins=np.linspace(-0.025,0.025, 50), normed=1, lw=2.5, color="purple", label="v2 Sim")
		if data2 is not None:
			ax.hist(data2["mean_psf_e1_sky"], histtype="step", bins=np.linspace(-0.025,0.025,50), normed=1, lw=2.5, color="steelblue", label="v02 Data")
		ax.set_xlabel("PSF Ellipticity $e^{PSF}_1$")

		ax2 = fig.add_subplot(132)
		plt.setp(ax2.get_yticklabels(), visible=False)
		ax2.set_xticks([-0.05,0,0.025,0.05])
		ax2.hist(data["mean_psf_e2_sky"], histtype="step", bins=np.linspace(-0.025,0.025, 50), normed=1, lw=2.5, color="purple", label="v2 Sim")
		if data2 is not None:
			ax2.hist(data2["mean_psf_e2_sky"], histtype="step", bins=np.linspace(-0.025,0.025,50), normed=1, lw=2.5, color="steelblue", label="v02 Data")
		ax2.set_xlabel("PSF Ellipticity $e^{PSF}_2$")

		ax3 = fig.add_subplot(133)
		plt.setp(ax2.get_yticklabels(), visible=False)
		plt.setp(ax3.get_yticklabels(), visible=False)
		ax3.set_xticks([0.01,0.02,0.03,0.04])
		e = np.sqrt(data["mean_psf_e1_sky"]*data["mean_psf_e1_sky"] + data["mean_psf_e2_sky"]*data["mean_psf_e2_sky"])
		ax3.hist(e, histtype="step", bins=np.linspace(0.,0.04, 50), normed=1, lw=2.5, color="purple", label="v2 Sim")
		if data2 is not None:
			e = np.sqrt(data2["mean_psf_e1_sky"]*data2["mean_psf_e1_sky"] + data2["mean_psf_e2_sky"]*data2["mean_psf_e2_sky"])
			ax3.hist(e, histtype="step", bins=np.linspace(0.,0.04,50), normed=1, lw=2.5, color="steelblue", label="v02 Data")
			matplotlib.rcParams['legend.fontsize']=10
			ax3.legend(loc="upper right")
		ax3.set_xlabel("PSF Ellipticity $|e^{PSF}|$")

		plt.subplots_adjust(hspace=0, wspace=0, left=0.01,right=0.97, top=0.62)
		plt.savefig("%s/psf_ellipticity-hist-v2sim-vs-y1v2data.png"%outdir)
		plt.close()

	if "psf_size" in names:
		print "-- PSF FWHM"
		plt.hist(data["mean_psf_fwhm"], histtype="step", bins=np.linspace(2.5,5.5, 50), normed=1, lw=2.5, color="purple", label="v2 Sim")
		if data2 is not None:
			plt.hist(data2["mean_psf_fwhm"], histtype="step", bins=np.linspace(2.5,5.5,50), normed=1, lw=2.5, color="steelblue", label="v02 Data")
			matplotlib.rcParams['legend.fontsize']=14
			plt.legend(loc="upper right")
		plt.xlabel("PSF Size $R^{PSF}_{FWHM}$")
		plt.xlim()

		plt.savefig("%s/psf_size-hist-v2sim-vs-y1v2data.png"%outdir)
		plt.close()

	if "size" in names:
		print "-- Radius"
		plt.hist(data["radius"], histtype="step", bins=np.linspace(0.,1.5, 50), normed=1, lw=2.5, color="purple", label="v2 Sim")
		if data2 is not None:
			plt.hist(data2["radius"], histtype="step", bins=np.linspace(0.,1.5,50), normed=1, lw=2.5, color="steelblue", label="v02 Data")
			matplotlib.rcParams['legend.fontsize']=14
			plt.legend(loc="upper right")
		plt.xlabel("Radius $R_{gpp}$")
		plt.xlim()

		plt.savefig("%s/radius-hist-v2sim-vs-y1v2data.png"%outdir)
		plt.close()








def histograms_vs_input(names, data, outdir="/home/samuroff/shear_pipeline/end-to-end/plots/v2.2/new/input", data2=None, thin=None):

	os.system("mkdir -p %s"%outdir)
	plt.close()

	print "Making plots:"

	if "flux" in names:
		print "-- flux"
		plt.hist(data["flux"], histtype="step", bins=np.linspace(0.0,4000, 50), normed=1, lw=2.5, color="purple", label="v2 Sim (input)")
		if data2 is not None:
			plt.hist(data2["mean_flux"], histtype="step", bins=np.linspace(0.0,4000,50), normed=1, lw=2.5, color="steelblue", label="v02 Data (measured)")
			plt.legend(loc="upper right")
		plt.xlabel("Galaxy Flux")
		plt.savefig("%s/flux-hist-v2sim_truth-vs-y1v2data.png"%outdir)
		plt.close()

	if "mag" in names:
		print "-- magnitude"
		plt.hist(data["mag"], histtype="step", bins=np.linspace(16,28, 50), normed=1, lw=2.5, color="purple", label="v2 Sim (input)")
		if data2 is not None:
			plt.hist(data2["mag_auto_r"], histtype="step", bins=np.linspace(16,28,50), normed=1, lw=2.5, color="steelblue", label="v02 Data (measured)")
			plt.legend(loc="upper left")
		plt.xlim(18,25)
		plt.xlabel("$r$-band Magnitude")
		plt.savefig("%s/mag-hist-v2sim_truth-vs-y1v2data.png"%outdir)
		plt.close()


	if "e" in names:
		print "-- Ellipticity"
		plt.subplots_adjust(wspace=0, hspace=0)
		#fig = plt.figure(1)
		ax=plt.subplot(131)
		plt.setp(ax.get_yticklabels(), visible=False)
		ax.set_xticks([-0.5,0,0.5,1])
		#ax.set_aspect(1)
		plt.hist(data["intrinsic_e1"], histtype="step", bins=np.linspace(-1,1, 50), normed=1, lw=2.5, color="purple", label="v2 Sim (input)")
		if data2 is not None:
			plt.hist(data2["e1"], histtype="step", bins=np.linspace(-1,1,50), normed=1, lw=2.5, color="steelblue", label="v02 Data (measured)")
		plt.xlabel("Ellipticity $e_1$")

		ax=plt.subplot(132)
		plt.setp(ax.get_yticklabels(), visible=False)
		ax.set_xticks([-0.5,0,0.5,1])
		#ax2.set_aspect(1)
		plt.hist(data["intrinsic_e2"], histtype="step", bins=np.linspace(-1,1, 50), normed=1, lw=2.5, color="purple", label="v2 Sim (input)")
		if data2 is not None:
			plt.hist(data2["e2"], histtype="step", bins=np.linspace(-1,1,50), normed=1, lw=2.5, color="steelblue", label="v02 Data (measured)")
		plt.xlabel("Ellipticity $e_2$")

		ax=plt.subplot(133)
		plt.setp(ax.get_yticklabels(), visible=False)
		ax.set_xticks([0.2,0.4,0.6,0.8,1])
		#ax3.set_aspect(1)intrinsic_
		e = np.sqrt(data["intrinsic_e1"]*data["intrinsic_e1"] + data["intrinsic_e2"]*data["intrinsic_e2"])
		plt.hist(e/(e.max()), histtype="step", bins=np.linspace(0.,1, 50), normed=1, lw=2.5, color="purple", label="v2 Sim (input)")
		if data2 is not None:
			e = np.sqrt(data2["e1"]*data2["e1"] + data2["e2"]*data2["e2"])
			plt.hist(e/e.max(), histtype="step", bins=np.linspace(0.,1,50), normed=1, lw=2.5, color="steelblue", label="v02 Data (measured)")
			matplotlib.rcParams['legend.fontsize']=10
			plt.legend(loc="upper right")
		plt.xlabel("Ellipticity $|e|$")

		plt.subplots_adjust(hspace=0, wspace=0, left=0.01,right=0.97, top=0.62)
		plt.savefig("%s/ellipticity-hist-v2sim_truth-vs-y1v2data.png"%outdir)
		plt.close()



	if "size" in names:
		print "-- Radius"
		plt.hist(data["hlr"], histtype="step", bins=np.linspace(0.,1.5, 50), normed=1, lw=2.5, color="purple", label="v2 Sim (input)")
		if data2 is not None:
			plt.hist(data2["radius"], histtype="step", bins=np.linspace(0.,1.5,50), normed=1, lw=2.5, color="steelblue", label="v02 Data (measured)")
			matplotlib.rcParams['legend.fontsize']=14
			plt.legend(loc="upper right")
		plt.xlabel("Radius $R_{gpp}$")
		plt.xlim(0.,1.5)

		plt.savefig("%s/radius-hist-v2sim_truth-vs-y1v2data.png"%outdir)
		plt.close()

	if "redshift" in names:
		print "-- redshift"
		plt.hist(data["cosmos_photoz"], histtype="step", bins=np.linspace(0.,1.6, 50), normed=1, lw=2.5, color="purple", label="v2 Sim (input)")
		if data2 is not None:
			plt.hist(data2["desdm_zp"], histtype="step", bins=np.linspace(0.,1.6,50), normed=1, lw=2.5, color="steelblue", label="v02 Data (measured)")
			matplotlib.rcParams['legend.fontsize']=14
			plt.legend(loc="upper right")
		plt.xlabel("Redshift")
		plt.xlim(0,1.6)

		plt.savefig("%s/redshift-hist-v2sim_truth-vs-y1v2data.png"%outdir)
		plt.close()
		











