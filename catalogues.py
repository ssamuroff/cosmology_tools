import numpy as np
import pyfits
import astropy.table as tb
import os.path as path

class catalogue:
    def __init__(self, catlist):

        tile = path.basename(catlist[0]).split("_")[0]
        if "DES" in tile:
            setattr(self,"tile",tile)
        
        self.ncat=len(catlist)
        self.bands=[]
        for i, catalogue_name in enumerate(catlist):
            catalogue = tb.Table.read(catalogue_name, format="fits")
            try:
                band = catalogue_name.split("_cat")[0].split("_")[-1]
                setattr(self,"%s"%band, catalogue)
                self.bands.append(band)
            except:
                setattr(self,"cat%d"%(i+1), catalogue)

            setattr(self, "nobj", len(catalogue))

            print "Loading %s, %d objects"%(catalogue_name, self.nobj)
        print "Loaded %d catalogues"%self.ncat

    def load_cosmos(self, loc="/share/des/disc2/image-sims/matched_cosmos_cat_25.2.fits"):
        cat = pyfits.getdata(loc)
        setattr(self,"cosmos_cat", cat)

    def plot_mag(self, band1, band2, newplot=0, colour="m.", xlim_lower=15.5, ylim_lower=15.5, xlim_upper=26.5, ylim_upper=26.5):
        """Generate a scatter plot of magnitudes in two bands for all objects in the tile."""
        import pylab as plt
        m1 = getattr(self,band1)["mag_auto"]
        m2 = getattr(self,band2)["mag_auto"]
        plt.plot(m1, m2, colour)
        plt.xlabel("$%s$"%band1)
        plt.ylabel("$%s$"%band2)
        plt.ylim(ylim_lower,ylim_upper)
        plt.xlim(xlim_lower, xlim_upper)

        plt.plot([0,30], [0,30], "k-")

    def assign_cosmos_ids(self, method="random", randomseed=90000):
        ncosmos = len(self.cosmos_cat)
        np.random.seed(randomseed)
        i = np.random.randint(0,ncosmos, self.nobj)

        ids= self.cosmos_cat[i]["ident"]

        newcol1 = tb.Column(i,"cosmos_index")
        newcol2 = tb.Column(ids, "cosmos_ident")

        for b in self.bands:
            exec "self.%s.add_column(newcol1)"%b
            exec "self.%s.add_column(newcol2)"%b

    def write_refcat(self,filename=None):

        if not filename:
            filename = "%s_cosmos_refcat.txt"%self.tile

        b=self.bands[0]
        cat = getattr(self, "%s"%b)

        x = cat["X_IMAGE"]
        y = cat["Y_IMAGE"]
        ind = cat["cosmos_index"]
        ident = cat["cosmos_ident"]

        out= np.vstack((x,y))
        out = np.vstack((out,ind))
        out = np.vstack((out,ident))

        np.savetxt(filename, out.T, header="x_image y_image cosmos_index cosmos_ident")
        print "Saving reference catalogue to %s"%filename









    def match(self, n1, n2):
         print "Selecting unique objects"
         cat1 = getattr(self, "cat%d"%n1)
         cat2 = getattr(self, "cat%d"%n2)
         if 'COADD_OBJECTS_ID' in cat1.columns:
             cat1.rename_column('COADD_OJECTS_ID', 'coadd_objects_id')
         if 'COADD_OJECTS_ID' in cat2.columns:
             cat2.rename_column('COADD_OJECTS_ID', 'coadd_objects_id')
         c1 = tb.unique(cat1, keys='coadd_objects_id')
         c2 = tb.unique(cat2, keys='coadd_objects_id')
         print "Sorting"
         c1 = c1.group_by("coadd_objects_id")
         c2 = c2.group_by("coadd_objects_id")
         
         tab = tb.join(c,d, join_type='inner')
    
         return tab
    
    def pixellise(self, ra, dec, observable, nside):
        import healpix_util as hu
        hpix = hu.HealPix('ring',nside)
        pixels = hpix.eq2pix(ra,dec)
        assert len(pixels)==len(observable)
        npix = hp.nside2npix(nside)
        notempty = np.unique(pixels)
        meanobs=np.zeros(npix)
        print "Averaging in %d (non empty) pixels"%len(notempty)
        for i in notempty:
            mask = (pixels==i)
            meanobs[i] = np.mean(observable[mask])
        return meanobs
