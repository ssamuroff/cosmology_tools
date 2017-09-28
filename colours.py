import numpy as np
import fitsio as fi
import pylab as plt
import scipy.interpolate as spi
plt.switch_backend("pdf")
import matplotlib 
matplotlib.rcParams["font.family"]="serif"

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

labels={
	"mean_z" : r"Mean Redshift $\bar{z}$",
	"snr" : r"Signal-to-Noise $S/N$"
}

lims={
	"mean_z" : [0.2,1.9]
}

class y1shear:
	def __init__(self,bpz=True,metacal=True,im3shape=True,gold=True, mof=True, verbose=True):

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
			if mof:
				self.bpz = fi.FITS("/global/cscratch1/sd/tvarga/bpz_matching/y1a1-gold-mof-badregion_BPZ.fits")[-1].read()
			else:
				self.bpz = fi.FITS("/global/cscratch1/sd/troxel/finaly1cats/mcal-y1a1-combined-griz-blind-v3-matched_BPZ.fits")[-1].read()


		self.verbose = verbose
		self.masks={}

	def get_from_cat(self,cat):

		attrs = dir(cat)
		for attr in attrs:
			# Copy all attributes of cat which are not functions
			is_function = hasattr( getattr(cat,attr), "__call__")
			if not is_function:
				print "Copying attribute %s"%attr
				val = getattr(cat, attr)
				setattr(self, attr, val)

	def get_f_R(self, nz_R, nz_B, split="colour", mask=[], xlim=[0.0,1.9], interpolate=True):
		select = self.get_type_split(split, mask=mask)

		if (len(mask)==0):
			mask = np.ones(self.bpz["coadd_objects_id"].size).astype(bool)

		red_frac = self.mcal["coadd_objects_id"][mask][select].size * 1.0 / self.mcal["coadd_objects_id"][mask].size
		print "Global red fraction : %3.2f"%red_frac

		# Load the nofz files from BPZ, and setup interpolators
		if self.verbose:
			print "Setting up interpolators: ",
		nR= np.loadtxt(nz_R).T
		binning_R=np.load(nz_R.replace("source.nz", "nz_source_zbin.npy")).T[1]
		nB = np.loadtxt(nz_B).T
		binning_B=np.load(nz_B.replace("source.nz", "nz_source_zbin.npy")).T[1]
		if self.verbose: print "red... ",
		pz_R = [ spi.interp1d(nR[0], n, kind="cubic") for n in nR[1:] ]
		if self.verbose: print "blue... ",
		pz_B = [ spi.interp1d(nB[0], n, kind="cubic") for n in nB[1:] ]
		nbins = len(nR[1:])

		zbins = np.linspace(xlim[0],xlim[1],100)
		if interpolate:
			z = (zbins[:-1]+zbins[1:])/2
		else:
			z = nR[0]
		yvec=[[] for n in nR[1:] ]

		for j in xrange(nbins):
			print "Processing bin %d/%d"%(j,nbins)
			NR = binning_R[binning_R==j].size
			NB = binning_B[binning_B==j].size
			prefactor = NR * 1.0 / (NR+NB)

			for i, z0 in enumerate(z):
				print i, z0,
				if interpolate:
					pr = pz_R[j](z0)
					pb = pz_B[j](z0)
				else:
					pr = nR[j+1][i]
					pb = nB[j+1][i]
				print pr, pb
				if (pr<0.5e-4):
					yvec[j].append( 0.0)
				else:
					yvec[j].append( prefactor * pr / (pr + pb) )

		return z, np.array(yvec) 


	def get_fofz(self, split="colour", mask=[], weights=[] ,xlim=[0.2,2.2]):
		if len(mask)==0:
			mask = np.ones(self.bpz.size).astype(bool)
		if len(weights)==0:
			weights = np.ones(self.bpz.size)

		bpz = self.bpz[mask]
		weights = weights[mask]
		isbulge = self.get_type_split(split, mask=mask)
		if len(isbulge)==0: return 0

		print "global bulge fraction : %3.3f"%(bpz["mean_z"][isbulge].size*1.0/bpz["mean_z"].size)

		zbins = np.linspace(xlim[0],xlim[1],90)
		z = (zbins[:-1]+zbins[1:])/2

		yvec=[]
		for i,(lower,upper) in enumerate(zip(zbins[:-1],zbins[1:])):
			print i, lower,upper, 
			select = (bpz["mean_z"]<upper) & (bpz["mean_z"]>lower)
			mean_wt =  weights[select].mean()
			mean_wtb =  weights[select][isbulge[select]].mean()

			# No of galaxies times the mean weight
			# If mean wt (bulges) == mean wt (all) the fraction reverts to N_bulge/N_all
			ngal = bpz["mean_z"][select].size 
			ngalb = bpz["mean_z"][select][ isbulge[select] ].size

			yvec.append([ngalb * mean_wtb, ngal  * mean_wt ])
			print ngalb*1.0/ngal

		return z, np.array(yvec)

	def get_type_split(self, split, mask=[]):
		"""Returns a selection mask for bukge galaxies"""
		if split in self.masks.keys(): return self.masks[split]

		if len(mask)==0:
			mask = np.ones(self.bpz.size).astype(bool)

		if (split=="colour"):
			x = self.gold["mof_flux_r"][mask]/self.gold["mof_flux_r"][mask].mean()
			colour = self.gold["mof_flux_g"][mask]/self.gold["mof_flux_z"][mask]

			param = [0.24,0.085]
			line = param[0]*x + param[1]
			select = (colour<line)
			
		elif (split=="bpz_template"):
			# Heymans dfn
			select =  (self.bpz["template_type"][mask]<1.0)
		elif (split=="im3shape_bord"):
			select = (self.i3s["bulge_fraction"][mask]==1)
		else:
			print "Unknown split type."
			print " Please choose from:"
			print "    - 'colour' for photometric colour based on MOF flux"
			print "    - 'bpz_template' for T_B"
			print "    - 'im3shape_bord' for im3shape bulge or disc fit"
			return []

		self.masks[split] = select
		return select

	def show_colourspace_division(self,mask, colours=["g","z"], save=True, colour_by=[None,"number",50], docolours=True, dolevels=["k","-",""], cmap="BuPu", pdf=True, title="", ylim=[0,0.7], xlim=[0,0.4], split_param=[0.24,0.085]):

		if len(mask)==0:
			print "Warning -- no masking (are you sure this is a good idea?)"
			mask = np.ones(self.gold.size).astype(bool)
		# Normalise the fluxes to the mean
		f0 = self.gold["mof_flux_r"][mask].mean()

		colour = self.gold["mof_flux_%s"%colours[0]][mask]/self.gold["mof_flux_%s"%colours[1]][mask]

		#Count the galaxies in grid cells
		print "Counting..."
		countsf,xbinsf,ybinsf = np.histogram2d(self.gold["mof_flux_r"][mask]/f0, colour, bins=500, range=[xlim,ylim], normed=1 )
		yf=(ybinsf[:-1]+ybinsf[1:])/2
		xf=(xbinsf[:-1]+xbinsf[1:])/2
		xxf,yyf=np.meshgrid(xf,yf)

		counts,xbins,ybins, = np.histogram2d(self.gold["mof_flux_r"][mask]/f0, colour, bins=80, range=[xlim,ylim], normed=1 )
		y=(ybins[:-1]+ybins[1:])/2
		x=(xbins[:-1]+xbins[1:])/2
		xx,yy=np.meshgrid(x,y)

		# Plot out whatever quantity is needed as a colour map in this space
		if docolours:
			print "Plotting colour map..."
			if (colour_by[1]=="number"):
			    #C = plt.contourf(xxf,yyf,countsf, cmap=cmap)
			    #plt.scatter(xxf.flatten(), yyf.flatten(), c=countsf.flatten(), cmap=cmap, marker=".", edgecolors="none")
				plt.hist2d(self.gold["mof_flux_r"][mask]/f0, colour, bins=260, range=[xlim,ylim], normed=1,cmap=cmap )
			else:
				print " ---- Colouring by %s"%colour_by[1]
				cat = getattr(self,colour_by[0])
				thin = colour_by[2]
				plt.scatter(self.gold["mof_flux_r"][mask][::thin]/f0, colour[::thin], c=cat[colour_by[1]][mask][::thin], edgecolors="none", cmap=cmap, marker=".")
				cb=plt.colorbar()
				plt.clim(lims[colour_by[1]])
				cb.set_label(label=labels[colour_by[1]],size=22)

		# Render the histograms as contours
		if isinstance(dolevels,list):
			print "Making histograms..."
			C = plt.contour(xx,yy,counts.T, 8, colors=dolevels[0],linestyles=dolevels[1], linewidth=.5,labels=dolevels[2])
			plt.legend(loc="upper left")
		plt.xlabel("$r-$band Flux $f_r$", fontsize=22)
		plt.ylabel("Photometric Colour $f_g/f_z$", fontsize=22)
		plt.xlim(xbins.min(),xbins.max())
		plt.ylim(ybins.min(),ybins.max())

		plt.title(title, fontsize=22)

		#green valley division
		xl = np.linspace(xlim[0],xlim[1],100)
		lin = split_param[0]*xl+split_param[1]
		plt.plot(xl,lin, color="forestgreen", ls="--", lw=2.5)

		if save:
			print "Saving..."
			extension = "pdf"*pdf + "png"*np.invert(pdf)
			filename = "/global/cscratch1/sd/sws/mof_colour_space_division-%s-colouredby_%s.%s"%(title,colour_by[1],extension)
			print filename
			plt.savefig(filename) 

def type_cut(pz):
    mask = np.zeros(pz.size)-9999
    mask[(pz['template_type']<1)]=1
    mask[(pz['template_type']>1)]=2
    return mask


def colour_cut(shapes, pz):
	# Hardcoded numbers
	a= [0.037,0.12,0.05,0.0]
	c0=[-0.1,-1.7,0.15,1.6]
	bins=[0.2,0.43,0.63,0.9,1.3]

	# Extract mags
	r = 30 - 2.5 * np.log10(shapes["flux_r"])
	z = 30 - 2.5 * np.log10(shapes["flux_z"])

	mask = np.zeros(shapes.size)-9999

	# Now loop over tomographic bins
	for b,(lower,upper) in enumerate(zip(bins[:-1],bins[1:])):
		select = (pz['mean_z']>lower) & (pz['mean_z']<upper)
		lin = (a[b] * r + c0[b])
		colour_select_r = ((r-z)>lin)
		mask[select & colour_select_r] = 1
		colour_select_b = ((r-z)<lin)
		mask[select & colour_select_b] = 2

	return mask

def colour_diagram(x0, y0, ls="-", colour="k", pdf=True, title="", ylim=[0,0.7], xlim=[0,0.4], split_param=[0.24,0.085]):
    #Count the galaxies in grid cells
    print "Counting..."
    counts,xbins,ybins, = np.histogram2d(x0, y0, bins=80, range=[xlim,ylim], normed=1 )
    y=(ybins[:-1]+ybins[1:])/2
    x=(xbins[:-1]+xbins[1:])/2
    xx,yy=np.meshgrid(x,y)
    print "Making histograms..."
    C = plt.contour(xx,yy,counts.T, 8, colors=colour,linestyles=ls, linewidth=.5)
    plt.xlabel("$r-$band Flux $f_r$", fontsize=22)
    plt.ylabel("$r-z$", fontsize=22)
    plt.xlim(xbins.min(),xbins.max())
    plt.ylim(ybins.min(),ybins.max())
    plt.title(title, fontsize=22)
    #green valley division
    if len(split_param)>0:
        xl = np.linspace(xlim[0],xlim[1],100)
        lin = split_param[0]*xl+split_param[1]
        plt.plot(xl,lin, color="forestgreen", ls="--", lw=2.5)
    return 0

def hist1d(C, pz, type_split=True, labels=True):
	N = float(C.size)
	if labels:
		lab = "All Galaxies"
	else:
		lab = None
	H, edges = np.histogram(C, bins=100, range=(-4,4))
	x = (edges[:-1]+edges[1:])/2
	plt.plot(x, H/N, color='purple', lw=2, label=lab)

	if type_split:
		#Define a split by BPZ template
		late = (pz['template_type']>1)
		early = (pz['template_type']<1)

		# Get the normalisations
		Ne = float(C[early].size)
		if labels:
			labe = "Early-type"
			labl = "Late-type"
		else:
			labe = None
			labl = None
		He,edgese = np.histogram(C[early], bins=100, range=(-4,4))
		xe = (edgese[:-1]+edgese[1:])/2

		Nl =float (C[late].size )
		Hl,edgesl = np.histogram(C[late], bins=100, range=(-4,4))
		xl = (edgesl[:-1]+edgesl[1:])/2

		plt.fill_between(xe, He/Ne, color='red', alpha=0.3, label=labe, linestyle=":")
		plt.fill_between(xl, Hl/Nl, color='royalblue', alpha=0.3, label=labl, linestyle="--")

	# Tweak the axes
	if labels:
		plt.legend(loc="upper left")
	plt.yticks(visible=False)
	print "Done"
	return None





