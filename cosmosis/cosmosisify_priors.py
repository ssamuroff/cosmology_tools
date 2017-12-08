from astropy.table import Table, Column
import numpy as np
import astropy.table
import sys
import glob
import fitsio
from scipy.interpolate import Rbf
import tools.emcee as mc

translation={
    "omegabh2"  : ("cosmological_parameters","omega_bh2"),
    "omegach2"  : ("cosmological_parameters","omega_mh2"),
    "theta"   : ("cosmological_parameters","100theta_mc"),
    "tau"  : ("cosmological_parameters","tau"),
    "omegak"  : ("cosmological_parameters","omega_k"),
    "mnu"  : ("cosmological_parameters","mnu"),
    "meffsterile"  : ("cosmological_parameters","meffsterile"),
    "w"  : ("cosmological_parameters","w"),
    "wa"  : ("cosmological_parameters","wa"),
    "nnu"  : ("cosmological_parameters","nnu"),
    "yhe"   : ("cosmological_parameters","yhe"),
    "logA"   : ("cosmological_parameters","ln(10^10as)"),
    "ns"  : ("cosmological_parameters","ns"),
    "r"  : ("cosmological_parameters","r") }

class chain(mc.chain):
	def __init__(self, filename, existing=None):
		if existing!=None:
			self.samples = existing.samples
			self.col_names = existing.col_names
			return None
		with open(filename) as fil:
			line = fil.readline().split()
		fil.close()
		names = [col.split("--")[-1] for col in line]
		print names

		self.samples = Table.read(filename, format="ascii", names=names)
		self.col_names = [col for col in names if ((col!="like") and (col!="weight"))]

	def setup_interpolator(self, cols=None):
		self.factors = {}
		for name in self.col_names: self.factors[name] = self.samples[name].max()
		if cols==None:
			cols = self.col_names

		args = [self.samples[name]/self.factors[name] for name in cols]
		args = np.vstack((args, self.samples["like"]))

		import pdb ; pdb.set_trace()

		Rbf()





class priors:
	def __init__(self, filename=None, cosmosis_format=False):
		self.range={}
		self.type={}
		self.section={}
		if filename!=None:
			if cosmosis_format:
				self.load_cosmosis_prior(filename)
			else:
				self.load(filename)

	def do_evidence_integral_kde(self, chain="/share/des/disc3/samuroff/cosmosis_format_chain-kids-450_fiducial-oldparam.txt", exclude=["omega_m", "sigma8_input"], max_iter=1e6):
		print "Setting up KDE interpolator."
		kde = NamedKDE(chain, exclude=exclude)
		total = 0
		num = 0
		converged = False
		check = []
		import pdb ; pdb.set_trace()
		for num in xrange(90000):
			#if num>=max_iter:
		#		converged = True

			names, values = self.draw_from_prior_space()
			like = kde(values, names)[0]
			print "iteration %d , like = %e"%(num,like)
			#total += like
			#num += 1
	
#			check.append(total/num)
#			if len(check)==100:
#				sigma = np.std(check)
#				print "Sigma = %f"%sigma
#				if sigma<tol:
#					converged = True
#				else:
#					check = check[50:]


	def load_cosmosis_prior(self, filename):
		section = None
		with open(filename) as fil:
			while True:
				line = fil.readline()
				if not line:
					break
				if line=="\n": 
					continue

				if "[" in line:
					section = line.split("[")[1].split("]")[0]
					continue
				elif "=" in line:
					name = line.split("=")[0].replace(" ", "")

				self.type[name] = line.split("=")[1].replace(" ", "")
				self.section[name] = section
				self.range[name] = [float(line.split()[-2]) ,float(line.split()[-1])]
				print name



	def load(filename, delta=False):
		
		with open(filename) as fil:
			while True:
				line = fil.readline()
				if not line:
					break
				if (line=="range"):
					continue
				name,lower,upper = line.split()
				if name not in translation.keys(): continue
				lower = float(lower)
				upper = float(upper)
				if (not delta) and (lower==upper):
					continue 
				new_name = translation[name][1]
				self.range[new_name] = (lower,upper)
				self.section[new_name] = translation[name][0]
				self.type[new_name] = "uniform"
				print "%s --> %s [%f, %f]"%(name, new_name, lower, upper)

	def centre_values(self):
		values = []
		names = []
		for name in self.range.keys():
			width = self.range[name][-1] - self.range[name][0]
			val =  self.range[name][0] + 0.5*width
			values.append(val)
			names.append(name)
		return names, values


	def draw_from_prior_space(self):
		values = []
		names = []
		for name in self.range.keys():
			names.append(name)
			irand = np.random.rand()

			width = self.range[name][-1] - self.range[name][0]
			val =  self.range[name][0] + irand*width

			values.append(val)
		return names, values


	def export(self, outfile):
		out = open(outfile, "wa")
		
		sections = np.unique([self.section[name] for name in self.range.keys()])
		for section in sections:
			lines="\n[%s]\n\n"%section
			for name in self.range.keys():
				if self.range[name][0]==self.range[name][1]: 
					continue
				line = "%s %s %f %f\n"%(name, self.type[name], self.range[name][0], self.range[name][1])
				sec = self.section[name]
				if sec!=section:
					continue
				else:
					lines+=line
			out.write(lines)

		out.close()

def nsample_multinest(filename):
	line = open(filename).readlines()[-3]
	if line.startswith('#nsample='):
		n = int(line[9:])
	else:
		n = None
	return n

class NamedKDE(object):
	def __init__(self, filename, exclude):
		from cosmosis.datablock import option_section, names
		from cosmosis.plotting.kde import KDE

		exclude += ['post', 'like']
		names = open(filename).readline().lstrip('#').split()
		names = [[i, name] for i, name in enumerate(names) if name.split("--")[-1] not in exclude]
		indices = np.array(names).T[0].astype(int)
		names =np.array(names).T[1].astype(str)
		n = nsample_multinest(filename)
		data = np.loadtxt(filename).T
		if n is not None:
			print "Using last {} rows of data".format(n)
			data = data[:,-n:]
		if "weight" in names:
			weight = data[np.argwhere(names=="weight")[0,0]]
			weight /= weight.max()
		else:
			weight = None

		indices = [index for index,name in zip(indices,names) if name!='weight' ]
		names = [name for name in names if name!='weight']
		self.data = data[indices]
		names = [name.split("--") for name in names]
		self.names = names

		self.kde = KDE(self.data, weights=weight)

	def reorder_values(self, names, values):
		check = np.array(self.names).T[1]
		if not (check==np.array(names)).all():
			#print "Will reorder parameters for interpolation."
			vals = []
			for name in check:
				#print name
				index = np.argwhere(np.array(names)==name)[0,0]
				vals.append(values[index])
			return vals
		else:
			return values

	def __call__(self, values, names):
		vals = self.reorder_values(names, values)
		return self.kde.normalize_and_evaluate(vals)
	

