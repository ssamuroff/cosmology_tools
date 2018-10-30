import numpy as np
import fitsio as fi
import os

class dvec:
	def __init__(self, filename):
		print('Reading from %s'%filename)
		self.exts = fi.FITS(filename)
		self.filename = filename

	def replace_dvec(self):

		# Make a copy of the input FITS file
		print(self.filename)
		os.system('cp %s %s'%(self.filename, self.filename.replace('.fits','-copy.fits')))
		outfits = fi.FITS('%s'%self.filename.replace('.fits','-copy.fits'), 'rw')

		# Work out how to slice the datavector
		edges = [0] + self.ends

		# Then cycle through the correlation types
		correlations = ['xip','xim','gammat','wtheta']
		for i, (lower,upper) in enumerate(zip(edges[:-1], edges[1:])):
			cname = correlations[i]

			newpts = self.pts[lower:upper]
			if len(newpts)==0:
				continue
			print('%s has %d datapoints.'%(cname, len(newpts)))

			outdat = self.exts[cname].read()
			outdat['VALUE'] = newpts
			outdat['BIN1'] = self.covorder[1][lower:upper]+1
			outdat['BIN2'] = self.covorder[2][lower:upper]+1
			outdat['ANGBIN'] = self.covorder[0][lower:upper]
			outfits[cname].write(outdat)

		outfits.close()
		print('Done.')

	def revert_corrs(self, names=['xip','xim','gammat','wtheta']):
		print('%d correlations to process.'%len(names))

		for name in names:
			outname = '0_%s.txt'%name

			d = self.exts[name].read()
			b1,b2 = d['BIN1'].astype(int)-1, d['BIN2'].astype(int)-1
			out = np.array([d['ANG'], b1, b2, d['VALUE']]).T
			np.savetxt(outname, out)

			print(outname)

		print('Done.')

	def reorder(self, covorder):
		self.covorder = covorder.T
		angbin = covorder.T[0]
		b1 = covorder.T[1]
		b2 = covorder.T[2]

		corrs = {'xip': (0, self.exts['xip'].read()),
		         'xim': (1, self.exts['xim'].read()),
		         'gammat': (2, self.exts['gammat'].read()),
		         'wtheta': (3, self.exts['wtheta'].read())
		         }

		corrnames = ['xip','xim','gammat','wtheta']

		hdr = self.exts['covmat'].read_header()
		self.ends = [hdr['STRT_%d'%j0] for j0 in [1,2,3]]
		self.ends.append(self.ends[-1]+len(self.exts['wtheta'].read()))
		i0,j0=0,0
		corr = 'xip'
		missing=[]
		cc = []

		for i,j,t in zip(b1,b2,angbin):
			if (i0==self.ends[j0]):
				j0+=1
				corr = corrnames[j0]
				print(corr)

			data = self.exts[corr].read()
			mask = (data['BIN1']==i+1) & (data['BIN2']==j+1) & (data['ANGBIN']==t)

			# If this is an autocorrelation then we can trivially swap the redshift bins			
			if (len(data["VALUE"][mask])==0) and (corr in ['xip','xim','wtheta']): 
				mask = (data['BIN1']==j+1) & (data['BIN2']==i+1) & (data['ANGBIN']==t)

			print(i,j,t,data["VALUE"][mask])

			cc.append(data["VALUE"][mask])

			if (len(data["VALUE"][mask])==0):
				missing.append((i,j,t,corr))

			i0+=1

		cc = np.concatenate(cc)
		self.pts = cc
		print('%d datapoints missing'%len(missing))




