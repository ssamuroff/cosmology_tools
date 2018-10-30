import argparse
import sys
import copy
import numpy as np

# Load the centroid values together
p0 = np.genfromtxt('%s/means.txt'%sys.argv[-1], dtype=[('name','S60'), ('value',float), ('std',float)])
l = np.genfromtxt('%s/lerr68.txt'%sys.argv[-1], dtype=[('name','S60'), ('value',float)])
u = np.genfromtxt('%s/uerr68.txt'%sys.argv[-1], dtype=[('name','S60'), ('value',float)])
# We're assuming the ordering is consistent here.
elow = abs(p0['value']-l['value'])
eupp = abs(p0['value']-u['value'])

nchain = p0[(p0['name']=='post')].size
ncols = len(p0)*1./nchain
nparam = len(u)*1./nchain

print('Found %d chain(s)'%nchain)
print('%d parameters per chain'%nparam)

# Split off the chains into separate arrays
P = np.array([ p0[i0*nparam:(i0+1)*ncols] for i0 in range(nchain) ])
El = np.array([ elow[i0*nparam:(i0+1)*nparam] for i0 in range(nchain) ])
Eu = np.array([ eupp[i0*nparam:(i0+1)*nparam] for i0 in range(nchain) ])


# Now just print everything.
# Sorry this is inefficient. But the arrays shouldn't be long enough for it to really matter.
for ichain in range(nchain):
	print('Chain %d'%(ichain+1))
	for name in P[ichain]['name']:
		if name=='post':
			continue
		mask = (P[ichain]['name']==name)
		print( '%s: %3.4f +%3.4f -%3.4f'%(name,P[ichain]['value'][mask],El[ichain][mask],Eu[ichain][mask]) )

print('Done')


