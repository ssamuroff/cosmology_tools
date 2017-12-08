from astropy.table import Table, Column
import numpy as np
import astropy.table
import sys
import glob

KEEP = "w  ns omegam  sigma8  H0  omegabh2 omeganuh2".split()

def add_extra_parameters(cat):
	"""
	Add additional DES systematics to the Planck chain. Since the data
	doesn't constrain these at all just draw from a flat distribution 
	over some allowed range.
	"""
	extraparam_file = args.extra
	extra = np.loadtxt(extraparam_file, dtype=[("param", "S50"), ("centre", "float64"), ("range", "float64")])
	print "Drawing random values for extra parameters ", extra["param"]

	nsamples = len(cat) 

	wt = cat["weight"]
	like = cat["like"]
	cat.remove_column("like")
	cat.remove_column("weight")

	for line in extra:
		param = line["param"]
		central_value = line["centre"]
		value_range = line["range"]

		samp = draw_from_range(nsamples, value_range, central_value)

		cat[param] = samp

	cat["weight"] = wt
	cat["like"] = like

	return cat

def draw_from_range(nsamples, value_range, central_value):
	vlower = central_value - value_range/2.0
	vupper = central_value + value_range/2.0
	return np.random.uniform(low=vlower, high=vupper, size=nsamples)

def write_cat(cat, filename, random=False, deweight=False):
	"""
	Write a catalog in the CosmoSIS format, optionally 
	converting a weight column into a repeated multiplicity 
	and reducing the weights to unity and also optionally 
	randomizing the catalog
	"""
	c = open(filename, 'w')
	c.write('#')
	c.write('    '.join(cat.colnames))
	c.write('\n')
	n = len(cat)
	order = np.arange(n, dtype=int)
	if random:
		np.random.shuffle(order)
	if deweight:
		weight = cat['weight']
		weight_col = cat.colnames.index("weight")
	for i in xrange(len(cat)):
		j = order[i]
		row = cat[j]
		if deweight:
			repeats = int(weight[j])
			row[weight_col] = 1.0
			for k in xrange(repeats):
				c.write('   '.join(str(s) for s in row))
				c.write('\n')
		else:
			c.write('   '.join(str(s) for s in row))
			c.write('\n')

	c.close()


def transform(cat):
	pairs = [
		('cosmological_parameters--w', 'w'),
		('cosmological_parameters--n_s', 'ns'),
		('cosmological_parameters--omega_m', 'omegam'),
		('cosmological_parameters--omega_b', 'omegabh2'),
		('cosmological_parameters--sigma8_input', 'sigma8'),
		('cosmological_parameters--s8', 'sigma8'),
		('cosmological_parameters--h0', 'H0'),
		('cosmological_parameters--omnuh2', 'omeganuh2'),
		]
	for new,old in pairs:
		if old in cat.colnames:
			cat.rename_column(old, new)
	KEEP = args.keep.split()
	if "H0" in KEEP:
		cat['cosmological_parameters--h0']/=100.0
	if "omegabh2" in KEEP:
		cat['cosmological_parameters--omega_b']/=cat['cosmological_parameters--h0']**2
	if "s8" in KEEP:
		cat['cosmological_parameters--s8'] = cat["cosmological_parameters--sigma8_input"] * np.sqrt(cat["cosmological_parameters--omega_m"]/0.3)

	#move weight and like to the end
	names = cat.colnames[2:] + cat.colnames[:2]
	cat = cat[names]
	return cat

def write_header(outfile, params):
	#put the like cols at the end
	line = "  ".join(params)
	outfile.write('#'+line+"\n")

def find_params(base):
	filename = base+".paramnames"
	names = ["weight", "like"]
	indices = [0,1]
	KEEP = args.keep.split()
	print "Keeping parameters %s"%KEEP 
	for i,line in enumerate(open(filename)):
		name = line.split()[0].rstrip('*')
		if name not in KEEP: continue
		if name=="s8": continue
		names.append(name)
		indices.append(i+2)
	return names,indices

def process_files(base, outfilename, randomize, deweight):
	filenames = glob.glob(base+"_[0-9].txt")
	params,indices = find_params(base)	
	chains = []
	for filename in filenames:
		print filename
		chains.append(np.loadtxt(filename, usecols=indices).T)
	chains = np.hstack(chains)
	cat=Table(rows=chains.T, names=params)
	cat = transform(cat)
	if args.extra:
		cat = add_extra_parameters(cat)
	write_cat(cat, outfilename, random=randomize, deweight=deweight)

def main():
	import argparse
	global args
	parser=argparse.ArgumentParser(description="Planck Format Chains -> Cosmosis Format")
	parser.add_argument("base", help="Root chain name including path")
	parser.add_argument("output", help="Output file name")
	parser.add_argument("--randomize", action='store_true', help="Randomize the ordering")
	parser.add_argument("--keep", default="w  ns omegam  sigma8  H0  omegabh2 omeganuh2 s8", help="Parameters to keep")
	parser.add_argument("--extra", default="", help="File with additional parameters to add.")
	parser.add_argument("--deweight", action='store_true', help="Convert to weight=1 chains by repeating lines (only for integer Planck weights i.e. not the _post_ ones)")
	if __name__ == '__main__':
		args = parser.parse_args()
		process_files(args.base, args.output, args.randomize, args.deweight)

	KEEP = args.keep.split()

main()


