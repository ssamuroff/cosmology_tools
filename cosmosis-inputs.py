import numpy as np
import scipy as sp

def generator:
	def nz(z,bins, filename):

		zedges = (z[1:]+z[:-1])/2.
		dz = zedges[1]-zedges[0]
		z = np.append(z[0]-dz,z)
		z = np.append(z, z[-1]+dz)

		even = np.arange(0,len(z), 2).astype(int)
		odd = np.arange(1,len(z), 2).astype(int)

		edges_lower = z[even]
		edges_upper = z[odd]

		out = np.array(edges_lower)
		out=np.vstack((out, edges_upper))

		for bin in bins:
			out = np.vstack((out, bin))

		np.savetxt(filename, out)
