import numpy as np
import fitsio as fi
# All of these catalogues should be pre-matched 

#im3shape
#i3s = fi.FITS("/global/cscratch1/sd/tvarga/im3shape_matching/y1a1-im3shape_v5_matched_v7.fits")[-1].read()
#select = (i3s["info_flag"]==0)
#
##metacal
#mcal = fi.FITS("y1a1-gold-mof-badregion.fits")[-1].read()
#select = (mc["flags_select"]==0)
#
##gold
#gold = fi.FITS("mcal-y1a1-combined-riz-blind-v4-matched.fits")[-1].read()
#
## photoz
#bpz = fi.FITS("/global/cscratch1/sd/tvarga/bpz_matching/y1a1-gold-mof-badregion_BPZ.fits")[-1].read()

class y1shear:
	def __init__(self,bpz=True,metacal=True,im3shape=True,gold=True):

		print "will load:"
		for name in ["bpz","metacal","im3shape","gold"]:
			print name, eval(name)
		#im3shape
		if im3shape:
			print "Loading im3shape catalogue"
			self.i3s = fi.FITS("/global/cscratch1/sd/tvarga/im3shape_matching/y1a1-im3shape_v5_matched_v7.fits")[-1].read()
			self.iselect = (self.i3s["flags_select"]==0)

		#metacal
		if metacal:
			print "Loading metacal catalogue"
			self.mcal = fi.FITS("/global/cscratch1/sd/sws/y1/shapes/mcal-y1a1-combined-riz-blind-v4-matched.fits")[-1].read()
			self.mselect = (self.mcal["flags_select"]==0)

		#gold
		if gold:
			print "loading gold catalogue"
			self.gold = fi.FITS("/global/cscratch1/sd/sws/y1/shapes/y1a1-gold-mof-badregion.fits")[-1].read()
		# photoz
		if bpz:
			print "loading photo-z catalogue"
			self.bpz = fi.FITS("/global/cscratch1/sd/tvarga/bpz_matching/y1a1-gold-mof-badregion_BPZ.fits")[-1].read()

	def get_from_cat(self,cat):

		attrs = dir(cat)
		for attr in attrs:
			# Copy all attributes of cat which are not functions
			is_function = hasattr( getattr(cat,attr), "__call__")
			if not is_function:
				print "Copying attribute %s"%attr
				val = getattr(cat, attr)
				setattr(self, attr, val)

	def get_fofz(self, split="colour", mask=[], weights=[] ):
		if len(mask)==0:
			mask = np.ones(self.bpz.size).astype(bool)
		if len(weights)==0:
			weights = np.ones(self.bpz.size)

		bpz = self.bpz[mask]
		weights = weights[mask]
		isbulge = self.get_type_split(split, mask=mask)

		print "global bulge fraction : %3.3f"%(bpz["mean_z"][isbulge].size*1.0/bpz["mean_z"].size)

		zbins = np.linspace(0.2,1.8,101)

		z = (zbins[:-1]+zbins[1:])/2

		yvec=[]
		for i,(lower,upper) in enumerate(zip(zbins[:-1],zbins[1:])):
			print i, lower,upper
			select = (bpz["mean_z"]<upper) & (bpz["mean_z"]>lower) 
			ngal = weights[select].sum()

			ngalb = weights[select & isbulge].sum()

			yvec.append([ngalb,ngal])

		return z, np.array(yvec)

	def get_type_split(self, split, mask=[]):
		"""Returns a selection mask for bukge galaxies"""

		if len(mask)==0:
			mask = np.ones(self.bpz.size).astype(bool)

		if (split=="colour"):
			x = self.gold["mof_flux_r"][mask]/self.gold["mof_flux_r"][mask].mean()
			colour = self.gold["mof_flux_g"][mask]/self.gold["mof_flux_z"][mask]

			line = -7*x + 0.7
			return (colour<line)
		elif (split=="bpz_template"):
			# Heymans dfn
			return (self.bpz["template_type"][mask]<2.0)
		elif (split=="im3shape_bord"):
			return (self.i3s["bulge_fraction"][mask]==1)
		else:
			print "unknown split type."
			print " Please choose from:"
			print "'colour' for photometric colour based on MOF flux"
			print "'bpz_template' for T_B"
			print "'im3shape_bord' for im3shape bulge or disc fit"

	def show_colourspace_division(self,mask, colours=["g","z"], coloured_by=[None,"number"], cmap="BuPu"):
		if (mask==None):
			mask = np.ones(self.gold.size).astype(bool)
		# Normalise the fluxes to the mean
		f0 = self.gold["mof_flux_r"][mask].mean()

		colour = self.gold["mof_flux_%s"%colours[0]][select]/self.gold["mof_flux_%s"%colours[1]][select]

		#Count the galaxies in grid cells
		countsf,xbinsf,ybinsf = np.histogram2d(self.gold["mof_flux_r"][mask]/f0,colour,bins=[np.linspace(0,0.4,500),np.linspace(0,0.7,500) ])
		yf=(ybinsf[:-1]+ybinsf[1:])/2
		xf=(xbinsf[:-1]+xbinsf[1:])/2
		xxf,yyf=np.meshgrid(xf,yf)

		counts,xbins,ybins, = np.histogram2d(self.gold["mof_flux_r"][mask]/f0,colour,bins=[np.linspace(0,0.4,80),np.linspace(0,0.7,80) ])
		y=(ybins[:-1]+ybins[1:])/2
		x=(xbins[:-1]+xbins[1:])/2
		xx,yy=np.meshgrid(x,y)

		# Render the histograms as contours
		plt.close()
		if (colour_by[1]=="number"):
			C = contourf(xxf,yyf,countsf, cmap=cmap)
		else:
			prin "Colouring by %s"%self,colour_by[1]
			cat = getattr(self,colour_by[0])
			plt.scatter(self.gold["mof_flux_r"][mask][::50]/f0, colour[::50], c=cat[colour_by[1]][::50], edgecolors="none")
			
		C=contour(xx,yy,counts, 8, colors='black', linewidth=.5)
		plt.xlabel("$r-$band Flux $f_r$")
		plt.ylabel("Photometric Colour $f_g/f_z$")
		plt.xlim(xbins.min(),xbins.max())
		plt.ylim(ybins.min(),ybins.max())

		#green valley division
		lin = -7*x+0.7
		plt.plot(x,lin, color="forestgreen", ls="--", lw=2.5)

		plt.savefig("/global/cscratch1/sd/sws/mof_colour_space_division.pdf") 