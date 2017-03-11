import numpy as np
import scipy as sp
import fitsio as fi
import glob
import tools.covariance as cv
import tools.n_of_z as nz

dirs={"xip":"shear_xi/", "xim":"shear_xi/"}
datavector_names={"xip":("theta","xiplus"), "xim":("theta","ximinus")}
realspace_lookup={"xip":True, "xim":True}
corr={"xip":("ee","xi+"), "xim":("ee","xi-")}
header_quantity_names = {"xip":("G+R","G+R"), "xim":("G-R","G-R"), "gammat":("G+R","GPR"), "wtheta":("GPR","GPR")}


class fits_object():
	def __init__(self,datavectors,covariance=None, nofz_shear=None, nofz_density=None, kernel_shear=None, kernel_density=None, verbosity=1, nongaussian=True):
		#Expects a list of datavectors eg ['xi+','xi-','gammat','wtheta']
		# And optionally a covariance wrapper object

		self.datavectors = datavectors
		self.update_kernel_names(kernel_shear,kernel_density)
		self.realspace = realspace_lookup[datavectors[0]]

		print "loading datavectors:"

		statistics=[]
		correlations=[]

		for name in datavectors:
			if verbosity==1:
				print name
			xname,basename = datavector_names[name]
			xpath = "%s/%s.txt"%(dirs[name],xname)
			x = np.loadtxt(xpath)
			setattr(self,xname,x)

			files = glob.glob("%s/%s*.txt"%(dirs[name],basename))
			for f in files:
				data = np.loadtxt(f)
				# Get the bin pairing from the filename
				i = int(f.split(".txt")[0].split("_")[-2])
				j = int(f.split(".txt")[0].split("_")[-1])

				setattr(self,"%s_%d_%d"%(name,i,j),data)
				if verbosity>1:
					print name, i, j


			c, st = corr[name]
			correlations.append(c)
			statistics.append(st)

		if (nofz_shear is not None):
			self.nofz_shear = nz.sv_pz(nofz_shear,"skynet")
			self.nofz_shear.load_from_txt()
			self.do_nz=True

		if (nofz_density is not None):
			self.nofz_density = nz.sv_pz(nofz_density,"skynet")
			self.nofz_density.load_from_txt()
			self.do_nz=True


		if covariance is not None:
			if verbosity>0:
				print "Loading covariance matrix"
			correlations = list(np.unique(correlations))
			# This function handles a lot of the fiddly details of the covariance matrix
			# You don't really want to delve into that
			# But the covariance wrapper should take care of all the symmetries, so self.cov should
			# have all elements, even where duplicated elsewhere in the matrix
			cov = covariance.extract_all(correlations,statistics=statistics, verbosity=verbosity, nongaussian=nongaussian)
			setattr(self,"cov",cov)
			setattr(self,"covariance_source", covariance)
			self.get_sampling()

	def get_sampling(self):
		"""Interpolate the datavectors to the same sampling in scale as the 
		covariance matrix."""
		print "Interpolating to match the covariance matrix"

		if self.covariance_source.realspace:
			k = 60.0 * 180.0/np.pi
		else:
			k=1.0
		
		for name in self.datavectors:
			xname,basename = datavector_names[name]
			files = glob.glob("%s/%s*.txt"%(dirs[name],basename))

			if self.covariance_source.realspace:
				xname="theta"
			else:
				xname="ell"
			for f in files:
				i = int(f.split(".txt")[0].split("_")[-2])
				j = int(f.split(".txt")[0].split("_")[-1])

				x = getattr(self.covariance_source,xname)
				x0 = getattr(self,xname) * k
				y0 = getattr(self,"%s_%d_%d"%(name,i,j))
				y = np.interp(x,x0,y0)
				setattr(self,"%s_%d_%d"%(name,i,j),y)

		setattr(self,xname,x)


	def export(self, filename, verbosity=1):
		fits = fi.FITS(filename, "rw", clobber=True)

		order={}
		dimension={}

		done_nz=[]

		if verbosity>0:
			print "Reformatting data for export to fits file."

		for obs in self.datavectors:
			if verbosity>0:
				print obs
			# Rearrange the data into the columns expected
			bin1, bin2, xbin, y, x = self.format_data(obs, verbosity)

			# Store this for the ordering of the covariance matrix
			order[obs] = zip(bin1,bin2)
			dimension[obs] = y.size

			if verbosity>1:
				print "writing to file..."
			# Extract the header information
			h = self.create_data_header(x,bin1,bin2,obs)
			out = self.write_for_export(bin1, bin2, xbin, y, x)
			fits.write(out, header=h)
			fits[-1].write_key("EXTNAME",obs)
			if verbosity>1:
				print "data"

			if self.do_nz:
				k1,k2 = self.header_kernel_names[obs]
				out1, out2 = self.construct_pz_hdus(obs)
		
				if k1 not in done_nz:
					fits.write(out1)
					fits[-1].write_key("EXTNAME",k1)
					done_nz.append(k1)

				if k2 not in done_nz:
					fits.write(out2)
					fits[-1].write_key("EXTNAME",k2)
					done_nz.append(k2)
	

		# Also include the covariance matrix if we have one loaded
		if hasattr(self, "cov"):
			if verbosity>1:
				print "covariance"
			cov = self.construct_covariance_matrix(order)
			h = self.create_covariance_header(dimension)

			fits.write(cov, header=h)
			fits[-1].write_key("EXTNAME","COV")

		if verbosity>1:
			print "done"

		fits.close()
		return 0

	def construct_pz_hdus(self, obs):
		if obs in ["xip","xim","cl_ee"]:
			pz1 = self.nofz_shear
			pz2 = self.nofz_shear
		elif obs in ["wtheta","cl_nn"]:
			pz1 = self.nofz_density
			pz2 = self.nofz_density
		elif obs in ["gammat","cl_ne"]:
			pz1 = self.nofz_shear
			pz2 = self.nofz_density

		dt1 = [ ('Z_LOW', '>f8'), ('Z_MID', '>f8'), ('Z_HIGH', '>f8')] + [("BIN%d"%(i+1), ">f8") for i in xrange(len(pz1.pz))]
		dt2 = [ ('Z_LOW', '>f8'), ('Z_MID', '>f8'), ('Z_HIGH', '>f8')] + [("BIN%d"%(i+1), ">f8") for i in xrange(len(pz2.pz))]

		out1 = np.zeros(pz1.z.size, dtype=dt1)
		out2 = np.zeros(pz2.z.size, dtype=dt2)

		if pz1.edges.size!=pz1.z.size+1:
			dz = pz1.z[1]-pz1.z[0]/2
			pz1.edges=np.linspace(pz1.z.min()-dz, pz1.z.max()+dz, pz1.z.size+1 )
		if pz2.edges.size!=pz2.z.size+1:
			dz = pz2.z[1]-pz2.z[0]/2
			pz2.edges=np.linspace(pz2.z.min()-dz, pz2.z.max()+dz, pz2.z.size + 1 )

		out1["Z_LOW"],out1["Z_HIGH"],out1["Z_MID"] = pz1.edges[:-1], pz1.edges[1:], pz1.z
		out2["Z_LOW"],out2["Z_HIGH"],out2["Z_MID"] = pz2.edges[:-1], pz2.edges[1:], pz2.z

		for i, p in enumerate(pz1.pz):
			out1["BIN%d"%(i+1)] = p

		for j, p in enumerate(pz2.pz):
			out2["BIN%d"%(j+1)] = p

		return out1,out2


	def construct_covariance_matrix(self, bin_pairs):
		"""Put together a covariance matrix already loaded into memory,
		   taking a row of blocks at a time"""

		for count1, obs in enumerate(self.datavectors):
			done=[]
			for count2, (i,j) in enumerate(bin_pairs[obs]):
				if [(i,j)] in done: continue
				row = self.get_row(bin_pairs, obs, i, j)

				if count1+count2==0:
					full_matrix = row
				else:
					full_matrix = np.vstack((full_matrix,row))

				done.append([(i,j)])

		return full_matrix



	def get_row(self, bin_pairs, obs1, i, j):
		
		for count1, obs2 in enumerate(self.datavectors):
			done=[]
			for count2, (k,l) in enumerate(bin_pairs[obs2]):
		
				if [(k,l)] not in done:
					correlation1, obs_name1 = corr[obs1]
					correlation2, obs_name2 = corr[obs2]
					c = tuple(list(correlation1)+list(correlation2))
					block = self.cov[(obs_name1,obs_name2)][c][(i,j,k,l)]

					if count1+count2==0:
						row= block
					else:
						row = np.hstack((row,block))
					done.append([(k,l)])

		return row


	def write_for_export(self, bin1, bin2, xbin, y, x):
		#There must be a neater way to do this.
		# We'll go with this for the moment because it works
		dt = np.dtype([('BIN1', '>i8'), ('BIN2', '>i8'), ('ANGBIN', '>i8'), ('VALUE', '>f8'), ('ANG', '>f8')])
		out=np.zeros(bin1.size, dtype=dt)
		out["BIN1"] = bin1
		out["BIN2"] = bin2
		out["ANGBIN"] = xbin
		out["VALUE"] = y
		out["ANG"] = x
		return out

	def create_covariance_header(self, dimension_lookup):
		header={}
		start=0
		for i,dat in enumerate(self.datavectors):
			header["NAME_%d"%i] = dat
			header["STRT_%d"%i] = start
			start += dimension_lookup[dat]

		return header


	def create_data_header(self,x,b1,b2,obs):
		units = "arcmin"
		# Assume arcmin for now
		nz1,nz2 = np.unique(b1).size, np.unique(b2).size
		return {"2PTDATA":True, "QUANT1": header_quantity_names[obs][0], "QUANT2": header_quantity_names[obs][1], "KERNEL_1": self.header_kernel_names[obs][0], "KERNEL_2": self.header_kernel_names[obs][1] ,"WINDOWS": "SAMPLE", "N_ZBIN_1": nz1, "N_ZBIN_2": nz2, "ANG_MIN": x.min(), "ANG_MAX": x.max(), "N_ANG": np.unique(x).size, "TUNIT5": units}

	def format_data(self,statistic,verbosity):
		bin1=[]
		bin2=[]
		angbin=[]
		value=[]
		ang=[]

		for datav_name in dir(self):
			if statistic not in datav_name:
				continue
			else:
				if verbosity>1:
					print datav_name

			i,j = int(datav_name.split("_")[-2]), int(datav_name.split("_")[-1])

			data = getattr(self,datav_name)
			nx = data.size

			xname = datavector_names[statistic][0]
			x = getattr(self,xname)

			bin1.append([i]*nx)
			bin2.append([j]*nx)
			angbin.append(list(np.linspace(0,nx-1,nx)))
			value.append(list(data))
			ang.append(list(x))

		return np.concatenate(bin1), np.concatenate(bin2), np.concatenate(angbin), np.concatenate(value), np.concatenate(ang)


	def update_kernel_names(self, kernel1, kernel2):
		self.header_kernel_names={}
		if kernel1 is None:
			kernel1 = "nz_shape"
		if kernel2 is None:
			kernel2 = "nz_pos"
		k={"e":kernel1,"n":kernel2}

		for obs_type in self.datavectors:
			c1,c2 = corr[obs_type][0]
			self.header_kernel_names[obs_type] = (k[c1],k[c2]) 

    		






















class generator:
	@staticmethod
	def nz(z,bins, filename):

		zedges = (z[1:]+z[:-1])/2.
		dz = zedges[1]-zedges[0]
		zedges = np.append(z[0] - dz/2,zedges)
		zedges = np.append(zedges, z[-1] + dz/2)

		edges_lower = zedges[:-1]
		edges_upper = zedges[1:]

		out = np.array(edges_lower)
		out=np.vstack((out, edges_upper))

		for bin in bins:
			out = np.vstack((out, bin))

		np.savetxt(filename, out.T)

def compare_cls(ell, cl1, cl2):
	import pylab as plt
	import numpy as np
	plt.subplot(211)
	plt.plot(ell, ell*ell*cl1)
	plt.plot(ell, ell*ell*cl2)
	plt.xlim(0.1,2500.)
	plt.xscale('log') ; plt.yscale('log')
	plt.subplot(212)
	diff = 1. - cl1/cl2
	i = np.argwhere(abs(diff) >= 0.01 )[0]
	ell0 = ell[i]
	print "1 percent deviation at ell=%e"%ell0
	i = np.argwhere(abs(diff) >= 0.05 )[0]
	ell1 = ell[i]
	print "5 percent deviation at ell=%e"%ell1
	i = np.argwhere(abs(diff) >= 0.15 )[0]
	ell2 = ell[i]
	print "15 percent deviation at ell=%e"%ell2
	i = np.argwhere(abs(diff) >= 0.25)[0]
	ell3 = ell[i]
	print "25 percent deviation at ell=%e"%ell3
	plt.plot(ell, diff)
	plt.axvline(ell0,linestyle='--', color='k')
	plt.axvline(ell1,linestyle='--', color='k')
	plt.axvline(ell2,linestyle='--', color='k')
	plt.axvline(ell3,linestyle='--', color='k')
	plt.axhline(0., color='k')
	plt.xlim(0.1,2500) ; plt.ylim(-0.1,0.1)
	plt.xscale('log')
	plt.subplot(211)
	plt.axvline(ell0,linestyle='--', color='k')
	plt.axvline(ell1,linestyle='--', color='k')
	plt.axvline(ell2,linestyle='--', color='k')
	plt.axvline(ell3,linestyle='--', color='k')
	plt.show()


class clres:
	def __init__(self, path=""):
		import glob, os
		ggl = glob.glob(os.path.join(path+"galaxy_position_shape_cross_cl_unbinned/bin*.txt"))
		shear = glob.glob(os.path.join(path+"galaxy_shape_cl_unbinned/bin*.txt"))
		pos = glob.glob(os.path.join(path+"galaxy_position_cl_unbinned/bin*.txt"))

		self.ell = np.loadtxt(os.path.join(path+"galaxy_position_shape_cross_cl_unbinned/ell.txt"))

		self.ggl={}
		self.shear={}
		self.pos={}

		for p in ggl:
			bin_name = os.path.basename(p).replace(".txt", "")
			self.ggl[bin_name] = np.loadtxt(p)
			print "Loaded GGL C(l) %s"%bin_name

		for p in pos:
			bin_name = os.path.basename(p).replace(".txt", "")
			self.pos[bin_name] = np.loadtxt(p)
			print "Loaded position-position C(l) %s"%bin_name

		for p in shear:
			bin_name = os.path.basename(p).replace(".txt", "")
			self.shear[bin_name] = np.loadtxt(p)
			print "Loaded shear-shear C(l) %s"%bin_name


	def compare_cls(self, cl2, bin1, bin2, cltype="shear"):
		import pylab as plt
		import numpy as np
		if (self.ell!=cl2.ell).all():
			print "ERROR: ell sampling does not match."

		sp1 = getattr(self, cltype)
		sp2 = getattr(cl2, cltype)
		bin = "bin_%d_%d"%(bin1,bin2)

		diff = 1. - sp1[bin]/sp2[bin]
		i = np.argwhere(abs(diff) >= 0.01 )[0]
		ell0 = self.ell[i]
		print "1 percent deviation at ell=%f"%ell0
		i = np.argwhere(abs(diff) >= 0.05 )[0]
		ell1 = self.ell[i]
		print "5 percent deviation at ell=%f"%ell1
		i = np.argwhere(abs(diff) >= 0.10 )[0]
		ell2 = self.ell[i]
		print "10 percent deviation at ell=%f"%ell2
		i = np.argwhere(abs(diff) >= 0.20)[0]
		ell3 = self.ell[i]
		print "20 percent deviation at ell=%f"%ell3



