import numpy as np
import glob, os, yaml, fitsio
import scipy.spatial as sps
import tools.diagnostics as di
import astropy.table as tb

reasons=["artifact", "smallstamp", "nogalaxy", "twogalaxies", "offcenter", "other"]
lookup={"artifact":1, "smallstamp":2, "nogalaxy":3, "twogalaxies":4, "offcenter":5, "other":6,"assumedbad":7}

levels={"xiping":[1,2,3,4,5,6], "trump":[1,2,3,5,6], "milverton":[1,2], "hewitt":[4],"sahlin":[5]}

class whistleblower:
    def __init__(self, validated="/home/samuroff/local/python/lib/python2.7/site-packages/tools/im3shape/cosmos/good_190117.tab", blacklist="/home/samuroff/local/python/lib/python2.7/site-packages/tools/im3shape/cosmos/bad_190117.tab"):
        print "Blacklisted COSMOS profiles from %s"%blacklist
        
        try:
            tab = tb.Table.read(blacklist, format="ascii", names=["ident", "reason"])
            self.reason, self.blacklisted_ids = tab["reason"], tab["ident"]
        except:
            tab = tb.Table.read(blacklist, format="ascii", names=["ident"])
            self.blacklisted_ids = tab["ident"]

        self.good_ids=np.loadtxt(validated)

    def research(self, catalogue):
        self.report = np.zeros(catalogue.res.size)

        for reason in reasons:
            blacklist = self.blacklisted_ids[self.reason==reason]
            mask = np.in1d(catalogue.truth["cosmos_ident"], blacklist)
            self.report[mask] = lookup[reason]

    def publish(self, level="xiping"):
        exclude = levels[level]
        exclude.append(7)
        return np.invert(np.in1d(self.report, exclude))

    def vindicate(self, catalogue):
        mask = np.in1d(catalogue.truth["cosmos_ident"], self.good_ids)
        self.report = np.ones(catalogue.res.size) * 7
        self.report[mask] = 0










def bad_tree(cat, bad_ids, truth=None, verbose=True, mask=None, ntiles=np.inf):
    if truth is None:
        truth = self.truth_path

    R=[]
    ids=[]

    if verbose:
        print "Setting up..."

    if mask is None:
        mask=np.ones_like(cat.res).astype(bool)

    cat.tiles = np.unique(cat.res["tilename"][mask]).astype("S12")
    object_tiles = cat.res["tilename"][mask].astype("S12")
    columns=['id', 'DES_id', 'cosmos_ident', 'cosmos_photoz', 'spectral_type', 'flags', 'star_flag', 'nearest_neighbour_pixel_dist', 'nearest_neighbour_index', 'true_g1', 'true_g2', 'intrinsic_e1', 'intrinsic_e2', 'intrinsic_e1_hsm', 'intrinsic_e2_hsm', 'hlr',  'ra', 'dec', 'mag', 'flux', 'nexp', 'mean_psf_e1', 'mean_psf_e2', 'mean_psf_fwhm']


    print "Searching for corrupted COSMOS profiles in %d truth tables."%cat.tiles.size
    for i, tile in enumerate(cat.tiles):
        if i > ntiles:
            break

        select = (object_tiles==tile)

        print i+1, tile, object_tiles[select].size

        filename = glob.glob("%s/%s*.fz"%(truth,tile))

        if len(filename)==0:
            print "Truth table is missing."
            continue

        if len(filename)>1:
            print "Warning - multiple truth tables found." 
            print "Will use %s"%filename[0]

        filename = filename[0]
        truth_table = fitsio.FITS(filename)[1].read(columns=columns)

        filename = glob.glob(("%s/%s*_cat.fits"%(truth,tile)).replace("truth", "OPS"))
        if len(filename)==0:
            print "SExtractor catalogue is missing."
            continue
        if len(filename)>1:
            print "Warning - multiple object catalogues tables found." 
            print "Will use %s"%filename[0]
        filename = filename[0]
        se_cat = fitsio.FITS(filename)[1].read()

        import pdb ; pdb.set_trace()
        select_tile = cat.res["tilename"]==tile

        # Frustratingly time consuming step to match up pixel positions with the correct object ids
        xy=np.vstack((se_cat["ALPHAWIN_J2000"], se_cat["DELTAWIN_J2000"]))
        xy0 = np.vstack((cat.truth["ra"][select_tile],cat.truth["dec"][select_tile]))
        tree = sps.KDTree(xy.T)
        results = tree.query(xy0.T, k=1)
        se_cat = se_cat[results[1]]

        # Set up a KD tree structure using only the bad profiles in this tile
        select_bad = np.in1d(truth_table["cosmos_ident"], bad_ids)
        xy_im_bad = np.vstack(( se_cat["XWIN_IMAGE"][select_bad], se_cat["YWIN_IMAGE"][select_bad] ))
        bad_tree = sps.KDTree(xy_im_bad.T)

        i3s_cuts = np.in1d(truth_table["DES_id"], cat.truth["DES_id"][select_tile])
        xy_im = np.vstack(( se_cat["XWIN_IMAGE"][i3s_cuts], se_cat["YWIN_IMAGE"][i3s_cuts] ))
        

        # Query it to find a match to each object
        R, ids = bad_tree.query(xy_im.T, k=2)
        import pdb ; pdb.set_trace()











        ra = truth_table["ra"][select_bad]
        dec = truth_table["dec"][select_bad]
        xy0 = np.vstack((ra[ra!=0.0], dec[dec!=0.0]))
        if verbose:
            print "Setting up KD tree (%d objects)"%ra.size,
        tree = sps.KDTree(xy0.T)

        # Query it to get a distance from every object in the final im3shape catalogue
        # to the nearest instance of a corrupted COSMOS profile
        if verbose:
            print "querying...",
        xy = np.vstack((cat.res["ra"][mask][select], cat.res["dec"][mask][select])) 
        result = tree.query(xy.T, k=1)[0] * 60 * 60 /0.27
        with open("corruption_dist.txt", "a") as f:
            out = np.vstack((cat.truth["DES_id"][mask][select], result))
            np.savetxt(f, out.T)
        R.append(result)
        ids.append(cat.res["coadd_objects_id"][mask][select])
        if verbose:
            print "Done %2.2f %2.2f"%(result.min(), result.max())

    b = np.zeros(np.concatenate(ids).size,dtype=[("coadd_objects_id", int),("dist", float)])
    b["coadd_objects_id"] = np.concatenate(ids)
    b["dist"] = np.concatenate(R)
    b,tmp = di.match_results(b,cat.res)

   
    return b["coadd_objects_id"], b["dist"]
