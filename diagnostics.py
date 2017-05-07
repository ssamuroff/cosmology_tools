import numpy as np
import astropy.table
#import pylab as plt
import glob, os 
import astropy.io.fits as pyfits
import pdb 
import fitsio as fi
import tools.arrays as arr
import scipy.optimize as optimise
from scipy.interpolate import Rbf

def text_to_fits(filename, keyword='fornax'):
    res = load_results0(keyword=keyword)
    fio.write(filename, res)
    print "Wrote merged results to %s"%filename

def rename(tilename):
    files = glob.glob("*fornax*.txt")
    for f in files:
        new = tilename+f
        print new
        os.system("mv %s %s"%(f, new))

def concatenate_results(outfile=None):
    if not outfile:
        print "Please specify an output filename."
    import glob
    files_main = glob.glob("*main*.txt")
    files_epoch = glob.glob("*epoch*.txt")

    catfile = outfile+"_main.txt"
    catfile_epoch = outfile+"_epoch.txt"
    
    f = open(files_main[0])
    cols = f.readline()
    add = f.readline()
    o = open(catfile, "wa")
    f.close()

    f = open(files_epoch[0])
    cols_epoch = f.readline()
    add_epoch = f.readline()
    oe = open(catfile_epoch, "wa")
    f.close()

    out=cols+" \n "+add+" \n "
    out_epoch=cols_epoch+" \n "+add_epoch+" \n "

    for i, fil in enumerate(files_main):
        f = open(fil)
        fe = open(files_epoch[i])
        tmp = f.readline()
        tmp = f.readline()
        data = f.read()
        tmp = fe.readline()
        tmp = fe.readline()
        data_epoch = fe.read()

        f.close()
        fe.close()

        out+=data
        out_epoch+=data_epoch
    #import pdb ; pdb.set_trace()

    o = open(catfile, "wa")
    o.write("%s"%out)
    o.close()

    o = open(catfile_epoch, "wa")
    o.write("%s"%out_epoch)
    o.close()

def load_epoch(epoch_path):
    print "Loading epoch catalogue from %s"%epoch_path
    import glob
    files=glob.glob("%s/*.fits"%epoch_path)
    e=[]                  
    for i, f in enumerate(files):
        e.append(pyfits.getdata(f))
        print i, f

    epoch = np.concatenate(np.array(e))

    print "Found %d epoch results"%epoch.size

    return epoch

def get_pixel_cols(meds_path):
    print "Loading pixel coordinate positions from %s"%meds_path
    import glob
    files=glob.glob("%s/DES*.fits*"%meds_path)

    x,y,cid,tile_identifier= [],[],[],[]        
    for i, f in enumerate(files):
        m = pyfits.open(f)["object_data"].data
        x.append(m["orig_col"].T[0])
        y.append(m["orig_row"].T[0])
        cid.append(m["id"])
        ti = np.array([i]*m["id"].size)
        tile_identifier.append(ti)
        print i, f

    tile_identifier = np.concatenate(tile_identifier)
    x = np.concatenate(np.array(x))
    y = np.concatenate(np.array(y))
    cid = np.concatenate(np.array(cid))

    print "Found %d epoch results"%cid.size

    out=np.zeros(cid.size,dtype=[("coadd_objects_id", int),("ix", float),("iy", float),("tile", int) ])
    out["ix"] = x
    out["iy"] = y
    out["coadd_objects_id"] = cid
    out["tile"] = tile_identifier

    return out






def load_results(res_path='None', keyword='fornax', format='txt', postprocessed=True, return_filelist=False, cols=None, ntot=-1, apply_infocuts=False, additional_cuts=[], match=[]):
    files=[]
    ind={}
    Nf=0
    i0=0
    if apply_infocuts:
        print "Applying info flag cuts on loading."
    if format=='txt':
        extra='.main'
    else:
        extra=''

    files = glob.glob('*%s*%s.%s'%(keyword,extra,format))
    if res_path!='None':
        files = glob.glob(os.path.join(res_path,'*%s*%s.%s'%(keyword,extra,format)))
    if len(files)==0:
        files = glob.glob('meds/*%s*%s.%s'%(keyword,extra,format))
    if len(files)==0:
        files = glob.glob('*%s*%s.%s'%(keyword,extra,format))

    if ntot!=-1:
        files=files[:ntot]

    if len(match)!=0:
        tiles = [os.path.basename(f)[:12] for f in files]
        select = np.array([t in match for t in tiles])
        files=np.array(files)[select]
        ntot = len(files)

    dt = fi.FITS(files[0])[1].read(columns=cols).dtype

  # Generate empty array to fill
  # Guess the size very roughly
    if not apply_infocuts:
        buff=50000
    else:
        buff=26000
    if ntot!=-1:
        res= np.zeros(ntot*buff, dtype=dt)
    else:
        res= np.zeros(len(files)*buff, dtype=dt)

    print "initial length of buffer : %d"%res.shape

            
    for f in files:
        tile = os.path.basename(f)[:12]
        if (len(match)>0):
            if  (tile not in match):
                #print "Excluding %s"%f
                #print "file exists but is not in the specified list of tiles."
                continue 
        if format.lower()=="fits":
            fits = fi.FITS(f)
            if cols:
                dat = fits[1].read(columns=cols)
            else:
                dat = fits[1].read()
            fits.close()
        elif format.lower()=="txt":
            if not cols:
                dat=np.array(astropy.table.Table.read(f, format="ascii")).astype(dt)
            else:
                dat = np.array(astropy.table.Table.read(f, format="ascii", include_columns=cols)).astype(dtype=dt)

        if apply_infocuts and 'info_flag' in res.dtype.names:
            dat = dat[dat["info_flag"]==0]

        if len(additional_cuts)>0:
            print "Applying additional cuts"
            count = additional_cuts.count("%s")
            cut = additional_cuts%tuple(["dat"]*count)
            exec cut
            dat = dat[cuts]

        nrows = len(dat)
        i1 = i0 + nrows
        if res[i0:i1].size<(i1-i0):
            print "Need more space - adding %d rows to the catalogue."%((i1-i0)*2)
            res = np.array(astropy.table.join(res, np.zeros((i1-i0)*2, dtype=dt), join_type="outer"))
        res[i0:i1] = dat
        ind[tile]=(i0,i1)
       # if res[i1]!=dat[-1]:
        #    print "ERROR: Data overflow. %d %d"%(i1, len(res))
        print Nf, tile, "(%d-%d)"%(i0,i1)
        Nf+=1
        i0 += nrows

    res = res[:i0]

    print 'Read results for %d objects from %d files.' %(len(res),len(files))

    if return_filelist:
        return res, files, ind
    else:
        return res

def load_results0(res_path='None', keyword='fornax'):
    files=[]
    files = glob.glob('*%s*.main.txt'%keyword)
    if res_path!='None':
        files = glob.glob(os.path.join(res_path,'*%s*.main.txt'%keyword))
    if len(files)==0:
        files = glob.glob('meds/*%s*.main.txt'%keyword)
    if len(files)==0:
        files = glob.glob('*%s*.main.txt'%keyword)
    res = np.array(astropy.table.Table.read(files[0], format='ascii')).astype(dtype=i3sdt)
    for f in files[1:]:
        res=np.concatenate((np.array(res).astype(i3sdt), np.array(astropy.table.Table.read(f, format='ascii')).astype(i3sdt)))
        print f

    print 'Read results for %d objects from %d files.' %(len(res),len(files))

    return res

#def load_truth(truth_path=None, keyword='DES'):
#    files=[]
#    if truth_path!=None:
#        files = glob.glob(os.path.join(truth_path,'*%s*-truth*.fits*'%keyword))
#        if len(files)==0: files = glob.glob(truth_path)
#    if len(files)==0:
#        files = glob.glob('../truth/*%s*-truth*.fits*'%keyword)
#    if len(files)==0:
#        files = glob.glob('*truth*/*%s*-truth*.fits*'%keyword)
#    if len(files)==0:
#        files = glob.glob('*%s*truth.*fits*'%keyword)
#    truth = pyfits.getdata(files[0],dtype=i3sdt)
#    
#    print 'Truth table contains %d objects' %len(truth)

#    return truth

def load_truth(truth_path=None, keyword='DES', match=None, apply_infocuts=True, cols=None, ind=None, res=None, faint=False, add_tilename_col=False):
    files=[]
    i0=0
    if truth_path!=None:
        files = glob.glob(os.path.join(truth_path,'*%s*-truth*.fits*'%keyword))
        if len(files)==0:
            files = glob.glob(truth_path)
    if len(files)==0:
        files = glob.glob('../truth/*%s*-truth*.fits*'%keyword)

    if match:
        list_res=[os.path.basename(m)[:12] for m in match]
        print "matching to tilelist ", list_res
            
        filelist=[]
        for f in files:
            if os.path.basename(f)[:12] in list_res:
                i=np.argwhere(np.array(list_res)==os.path.basename(f)[:12])[0,0]
                filelist.append(f)

        filelist = np.array(filelist)

  #      assert filelist.size == len(list_res)

        #filelist = np.array([np.array([ f, np.argwhere(np.array(list_res)==os.path.basename(f)[:12])[0,0] ]) for f in files if os.path.basename(f)[:12] in list_res])
    else:
        filelist=files

    if not faint:
        extension=1
        bookmark="DES_id"
    else:
        extension="subdetection_objects"
        bookmark="cosmos_id"


    dt = fi.FITS(filelist[0])[extension].read(columns=cols).dtype

    if apply_infocuts:
        buff=26000
    else:
        buff = 50000
    if faint: buff*=1000
    truth = np.zeros(len(filelist)*buff, dtype=dt)
    if add_tilename_col and ("tilename" not in dt.names):
        truth = arr.add_col(truth, "tilename", len(filelist)*buff*["DES0000+0000"])
        truth = arr.add_col(truth, "tile", len(filelist)*buff*[-9999])
    for i, f in enumerate(filelist):
        tile = os.path.basename(f)[:12]
        fits = fi.FITS(f)
        if cols:
            dat = fits[extension].read(columns=cols)
        else:
            dat = fits[extension].read()
        fits.close()

        if ind!=None and res!=None:
            mask = np.in1d(dat["DES_id"], res["coadd_objects_id"][ind[tile][0]:ind[tile][1]])
            dat = dat[mask]


        if add_tilename_col and ("tilename" not in dat.dtype.names):
            dat = arr.add_col(dat, "tilename", dat.shape[0]*[tile])
            dat = arr.add_col(dat, "tile", dat.shape[0]*[i])

        nrows = len(dat)
        i1 = i0 + nrows

        while truth[i0:i1].size!=(i1-i0):
            print "Need more space - adding %d rows to the truth table."%((i1-i0)*2)
            truth = np.array(astropy.table.join(truth, np.zeros((i1-i0)*2, dtype=dt), join_type="outer"))
        truth[i0:i1] = dat
        print i+1, tile, "%d/%d"%(i1, len(truth))
        i0 += nrows 
    #truth = np.concatenate((np.array(truth), np.array(astropy.table.Table.read(f, format="fits"))))
    
                
    truth = truth[truth[bookmark]!=0]

    print 'Truth table contains %d objects' %len(truth)

    return truth

def get_star_mask(tr):
    return tr['star_flag']!=1

def info_cuts(res):
    return res[res['info_flag']==0]

def get_bord_results(disc,bulge):
    bm=bulge['likelihood']>disc['likelihood']
    dm=bulge['likelihood']<disc['likelihood']
    if np.all(dm): print 'Warning: The disc fit is always chosen for these results.'
    elif np.all(bm): print 'Warning: The bulge fit is always chosen for these results.'
    
    # Apply selection
    res_bulge = bulge[bm]
    res_disc = disc[dm]
    # Combine
    res_bord = np.concatenate((res_disc, res_bulge))
    # Resort
    order = np.argsort(res_bord['id'])
    return res_bord[order]

def match_results(res,tr, table3=None, name1="coadd_objects_id", name2="coadd_objects_id", name3=None, unique=False, verbose=True):

    if table3 is not None:
        third_dataset=True
        if verbose:
            print "Found 3 catalogues to match"
    else:
        third_dataset=False

    if unique:
        if verbose:
            print "Enforcing unique entries."
        un,ind=np.unique(res[name2], return_index=True)
        res=res[ind]
        un,ind=np.unique(tr[name1], return_index=True)
        tr=tr[ind]
        if third_dataset:
            un,ind = np.unique(table3[name3], return_index=True)
            table3 = table3[ind]

    if third_dataset:
        table3 = table3[np.argsort(table3[name3])]

    if verbose:
        print "Sorting..."
    tr = tr[np.argsort(tr[name1])]
    res = res[np.argsort(res[name2])]

    if (res.shape==tr.shape) and (np.all(res[name2]==tr[name1])):
        print 'Matched %d objects'%len(tr[name1])
    else:
        if verbose:
            print "Matching..."

        tr = tr[np.in1d(tr[name1], res[name2])]
        res = res[np.in1d(res[name2], tr[name1])]

        if third_dataset:
            if verbose:
                print "Matching to third array"
            table3 = table3[np.in1d(table3[name3], res[name2])]
            table3 = table3[np.in1d(table3[name3], tr[name1])]

            tr = tr[np.in1d(tr[name1], table3[name3])]
            res = res[np.in1d(res[name2], table3[name3])]

            return res, tr, table3

    return res,tr 

def fast_match_results(table1, table2, name1="coadd_objects_id", name2="coadd_objects_id"):
    print "Sorting"
    indices2 = np.argsort(table2[name2])
    indices1 = np.argsort(table1[name1])

    print "Matching"
    indices = np.searchsorted(table1[name1],table2[name2])

    return table1[indices1], table2[indices2][indices]


def truth_histograms(tr, mask, mode='show'):
    dat = tr[mask]
    g1,g2 = dat['e1'], dat['e2']

    plt.hist(g1, bins=np.linspace(-0.1,0.1,50), histtype='step', color='m', label='$g_1$')
    plt.hist(g2, bins=np.linspace(-0.1,0.1,50), histtype='step', color='b', label='$g_2$')
    plt.xlim(-0.11,0.11)
    plt.axvline(0.0,color='k')
    plt.legend(loc='upper right')
    plt.xlabel('$g_i^{true}$')

    plt.hist(dat['intrinsic_e1'][dat['intrinsic_e1']!=0.0], bins=np.linspace(-1,1,50), histtype='step', color='m', label='$e_1$')
    plt.hist(dat['intrinsic_e2'][dat['intrinsic_e2']!=0.0], bins=np.linspace(-1,1,50), histtype='step', color='b', label='$e_2$')
    plt.xlim(-1,1)
    plt.axvline(0.0,color='k')
    plt.legend(loc='upper right')
    plt.xlabel('$e_i^{int}$')

    star_mask= (dat['star_flag']!=1)
    plt.hist(dat['flux'][star_mask], bins=np.linspace(0.,2000, 50), histtype='step', color='m', label='galaxies')
    if not np.all(star_mask):
        plt.hist(dat['flux'][np.invert(star_mask)], bins=np.linspace(0.,2000, 50), histtype='step', color='b', label='stars')
        plt.legend(loc='upper right')
    plt.xlabel('flux')

    star_mask= (dat['star_flag']!=1)
    plt.hist(dat['hlr'][star_mask], bins=np.linspace(0.,2.5, 50), histtype='step', color='m')
    plt.xlabel('hlr')

    plt.hist(dat['mean_psf_fwhm'], bins=np.linspace(0,0.75, 50), histtype='step', color='m')
    plt.xlabel('mean PSF FWHM')

    plt.hist(dat['mean_psf_e1_pixel'], bins=np.linspace(-0.1,0.1, 50), histtype='step', color='m', label='$e_1$')
    plt.hist(dat['mean_psf_e2_pixel'], bins=np.linspace(-0.1,0.1, 50), histtype='step', color='b', label='$e_2$')
    plt.legend(loc='upper right')
    plt.xlabel('mean $e_{i}^{psf}$')
    plt.xlim(-0.05,0.05)
    plt.axvline(0., color='k')

    plt.hist(dat['cosmos_photoz'][star_mask], bins=np.linspace(0.,2.0, 50), histtype='step', color='m')
    plt.xlabel('true redshift')

def binned_plots(res, savedir=None):
    print "plotting mean e vs"
    cols = res.dtype.names
    for c in cols:
        if not isinstance(res[c][0], float) and (c!="stamp_size"):
            continue
        if "covmat" in c or (len(np.unique(res[c]))<3.):
            continue
        print c
        
        meane_binned_plot(res, binned_by=c, savedir=savedir)
        

    print "done"



def meane_binned_plot(res, mask=None, binned_by="snr", savedir=None):
    lims={"snr":(10,150), "mean_rgpp_rp":(1.0,3.0), "round_snr": (10,300), "bulge_a":(400,5000), "disc_a":(400,5000), "mean_psf_fwhm":(2.9,5.), "mean_psf_e1": (-0.03, 0.03), "mean_psf_e2":(-0.03,0.03), "levmar_iterations":(0.,120), "mean_flux":(0,10000), "mean_mask_fraction":(0.,1.), "chi2_pixel": (0.8,1.3), "disc_flux": (0,10), "bulge_flux": (0.,10), "stamp_size":(20,300), "e1": (-1,1), "e2":(-1,1), "ra_as":(-2.,2.), "dec_as" : (-2,2), "radius":(0.,2.5)}

    if not mask:
        mask = np.ones_like(res['e1']).astype(bool)
    e1 = res[mask]['e1']
    e2 = res[mask]['e2']
    q = res[binned_by]
    sel = (e1>-1.) & (e1<1) & (e2>-1.) & (e2<1) & (np.isfinite(q))
    e1 = e1[sel]
    e2 = e2[sel]
    q = q[sel]
    if binned_by in lims.keys():
        bins = np.linspace(lims[binned_by][0], lims[binned_by][1], 9)
    else:
        bins = np.linspace(res[binned_by].min(), res[binned_by].max(), 9)

    print "mean: %f (%f, %f)" %(q.mean(), q.min(), q.max()) 

        
    meane1=[]
    meane2=[]
    vare1=[]
    vare2=[]
    meanq=[]
    varq=[]
    
    for i in xrange(len(bins)-1):
        sel1 = (q>bins[i]) & (q<bins[i+1])
        if np.isfinite(np.mean(e1[sel1])) and np.isfinite(np.mean(e2[sel1])) and np.isfinite(np.mean(q[sel1])):
            meane1 += [np.mean(e1[sel1])]
            meane2 += [np.mean(e2[sel1])]
            vare1 += [np.std(e1[sel1])/np.sqrt(len(e1[sel1]))]
            vare2 += [np.std(e2[sel1])/np.sqrt(len(e2[sel1]))]
            meanq.append(np.mean(q[sel1]))
            varq.append(np.std(e1[sel1])/np.sqrt(len(e1[sel1])))

    print bins
    print meane1
    print meane2
    if len(meane1)>0:
        plt.errorbar(meanq, meane1, yerr=vare1, fmt="o", color="m", label="$e_1$")
        plt.errorbar(meanq, meane1, yerr=vare1, fmt="-", color="m")
        plt.errorbar(meanq, meane2 , yerr=vare2, fmt="o", color="b", label="$e_2$")
        plt.errorbar(meanq, meane2, yerr=vare2, fmt="-", color="b")
        plt.legend()
        plt.xlabel(binned_by)
        plt.axhline(0.0,color='k')
        plt.ylabel("mean $e_i$")
    
        if savedir==None:
            plt.show()
        else:
            plt.savefig(savedir+"meane_binned_%s"%binned_by)
            plt.close()

def stabilised_fit(x,y,wt):
    print "Stabilising fit to get finite covariance (perhaps consider using more datapoints?)"
    import sys
    x_new = list(x) + list(x)[-1:]
    y_new = list(y) + list(y)[-1:]
    wt_new = list(wt) + [sys.float_info.epsilon]

    return np.polyfit(x_new,y_new, 1, w=wt_new, cov=True, full=False)




def truth_plots(res, tr, mask=None, nbins=5, mode='show', savedir=None, true_shape=False, psf_leakage=True, bin='', return_vals=False, use_intrinsic=False, equal_number_bins=True):
    if not mask:
        mask = np.ones_like(tr['true_g1']).astype(bool)
    g1 = tr[mask]['true_g1']
    g2 = tr[mask]['true_g2']
    sel = (g1>-1.) & (g1<1) & (g2>-1.) & (g2<1)
    g1=g1[sel]
    g2=g2[sel]
    if use_intrinsic:
        e1 = tr[mask]["intrinsic_e1"]
        e2 = tr[mask]["intrinsic_e2"]
    else:
        elab="e"
        e1 = res[mask]['%s1'%elab][sel]
        e2 = res[mask]['%s2'%elab][sel]
    e1int = tr[mask]['intrinsic_e1'][sel]
    e2int = tr[mask]['intrinsic_e2'][sel]
    e1psf = res[mask]['mean_psf_e1_sky']
    e2psf = res[mask]['mean_psf_e2_sky']

    if not equal_number_bins:
        bins = np.linspace(-0.03, 0.03, nbins+1)
    else:
        bins = find_bin_edges(g1, nbins)
        
    meane1=[] ;meane2=[]
    vare1=[] ;vare2=[]
    meang1=[] ;meang2=[]
    varg1=[] ; varg2=[]

    meane1_cross=[]
    meane2_cross=[]
    vare1_cross=[]
    vare2_cross=[]
    
    for i in xrange(len(bins)-1):
        sel1 = (g1>bins[i]) & (g1<bins[i+1])
        sel2 = (g2>bins[i]) & (g2<bins[i+1])
        meane1 += [np.mean(e1[sel1])]
        meane2 += [np.mean(e2[sel2])]
        meane1_cross += [np.mean(e1[sel2])]
        meane2_cross += [np.mean(e2[sel1])]
        vare1 += [np.std(e1[sel1])/np.sqrt(len(e1[sel1]))]
        vare2 += [np.std(e2[sel2])/np.sqrt(len(e2[sel2]))]
        vare1_cross += [np.std(e1[sel2])/np.sqrt(len(e1[sel2]))]
        vare2_cross += [np.std(e2[sel1])/np.sqrt(len(e2[sel1]))]
        meang1 += [np.mean(g1[sel1])]
        meang2 += [np.mean(g2[sel2])]
        varg1 += [np.std(g1[sel1])/np.sqrt(len(g1[sel1]))]
        varg2 += [np.std(g2[sel2])/np.sqrt(len(g2[sel2]))]
     #   plt.axvline(bins[i],color='k',linestyle = '--')
    
    meang1 = np.array(meang1)
    meang2 = np.array(meang2)
    meane1 = np.array(meane1)
    meane2 = np.array(meane2)
    varg1 = np.array(varg1)
    varg2 = np.array(varg2)
    vare1 = np.array(vare1)
    vare2 = np.array(vare2)
    de1 = meane1-meang1
    de2 = meane2-meang2
    vard1 = np.sqrt(vare1*vare1 + varg1*varg1)
    vard2 = np.sqrt(vare2*vare2 + varg2*varg2)
    
    sel = (g1>-1.) & (g1<1) & (g2>-1.) & (g2<1) 
    p11, covm11 = np.polyfit(meang1,de1, 1, w=1./vard1/vard1, cov=True, full=False)

    if not np.isfinite(covm11).all():
        p11,covm11 = stabilised_fit(meang1,de1,1.0/vard1/vard1)
    m11,c11 = p11
    p22, covm22 = np.polyfit(meang2,de2, 1, w=1./vard2/vard2, cov=True, full=False)
    if not np.isfinite(covm22).all():
        p22,covm22 = stabilised_fit(meang2,de2,1.0/vard2/vard2)
    m22,c22 = p22
#    m11 -=1.
#    m22 -=1.
    
    p12, covm12 = np.polyfit(meang1,meane2_cross, 1, w=1./vare1/vare1, cov=True, full=False)
    if not np.isfinite(covm12).all():
        p12,covm12 = stabilised_fit(meang1,meane2_cross,1.0/vare1/vare1)
    m12,c12 = p12
    p21, covm21 = np.polyfit(meang2,meane1_cross, 1, w=1./vare2/vare2, cov=True, full=False)
    if not np.isfinite(covm21).all():
        p21,covm21 = stabilised_fit(meang2,meane1_cross,1.0/vare2/vare2)
    m21,c21 = p21

    print 'm11=%f +- %f'%(m11,covm11[0,0]**0.5)
    print 'm22=%f +- %f'%(m22,covm22[0,0]**0.5)
    print 'm12=%f +- %f'%(m12,covm12[0,0]**0.5)
    print 'm21=%f +- %f'%(m21,covm21[0,0]**0.5)
    print 'c11=%f +- %f'%(c11,covm11[1,1]**0.5)
    print 'c22=%f +- %f'%(c22,covm22[1,1]**0.5)
    print 'c12=%f +- %f'%(c12,covm12[1,1]**0.5)
    print 'c21=%f +- %f'%(c21,covm21[1,1]**0.5)

    print "Total number of galaxies : %d"%len(res)

    if "is_bulge" in res.dtype.names:
        fb = 1.0*sum(res['is_bulge'])/len(res)
        print "Bulge fraction : %f" %(fb)
    else:
        fb =0.0

    if mode=="show":
        plt.errorbar(meang1,de1,yerr=vare1, color='m', fmt = 'o', label='$i,j=1,1$')
        plt.errorbar(meang2,de2,yerr=vare2, color='b', fmt='o', label='$i,j=2,2$')
        plt.plot(meang1,meang1*m11+c11, color='m')
        plt.plot(meang2,meang2*m22+c22, color='b')
    
        plt.axvline(0.,color='k')
        plt.axhline(0.,color='k')
        plt.legend(loc='upper right')
        plt.xlabel('$<g_i>$')
        plt.ylabel('$<e_j>-<g_j>$')
        plt.xlim(bins[0],bins[-1])
    if savedir:
        f = os.path.join(savedir, "%s_e-vs-g-truth-plot-auto.png"%bin)
        plt.savefig(f)
        plt.close()
    if mode=='show':
        plt.show()
        plt.close()

    if mode=="show":
        plt.errorbar(meang1,meane2_cross,yerr=vare2_cross, color='m', fmt = 'o', label='$i,j=1,2$')
        plt.errorbar(meang2,meane1_cross,yerr=vare1_cross, color='b', fmt='o', label='$i,j=2,1$')
        plt.plot(meang1,meang1*m12+c12, color='m')
        plt.plot(meang2,meang2*m21+c21, color='b')
        
        plt.axvline(0.,color='k')
        plt.axhline(0.,color='k')
        plt.legend(loc='upper right')
        plt.xlabel('$<g_i>$')
        plt.ylabel('$<e_j>$')
        plt.xlim(bins[0],bins[-1])
    if savedir:
        f = os.path.join(savedir, "%s_e-vs-g-truth-plot-cross.png"%bin)
        plt.savefig(f)
        plt.close()
    if mode=='show':
        plt.show()
        plt.close()


    if not isinstance(return_vals,bool):
        if return_vals=="m":
            return (m11+m22)/2.0, np.sqrt(covm11[0,0]+covm22[0,0])
        elif return_vals=="m1":
            return m11, np.sqrt(covm11[0,0])
        elif return_vals=="m2":
            return m22, np.sqrt(covm22[0,0])
        elif return_vals=="c":
            return (c11+c22)/2.0, np.sqrt(covm11[1,1]+covm22[1,1])
        elif return_vals=="c1":
            return c11, np.sqrt(covm11[1,1])
        elif return_vals=="c2":
            return c22, np.sqrt(covm22[1,1])

    if true_shape:
        meane1,meane2=[],[]
        meang1,meang2=[],[]
        vare1,vare2=[],[]
        varg1,varg2=[],[]
        meane1int,meane2int=[],[]
        vare1int,vare2int=[],[]
        true_sheared_shape=[[],[]]

        bins = np.linspace(-0.2, 0.2, nbins+1)
        for i in xrange(len(bins)-1):
            sel1 = (e1int>bins[i]) & (e1int<bins[i+1])
            sel2 = (e2int>bins[i]) & (e2int<bins[i+1])
            meane1 += [np.mean(e1[sel1])]
            meane2 += [np.mean(e2[sel2])]
            vare1 += [np.std(e1[sel1])/np.sqrt(len(e1[sel1]))]
            vare2 += [np.std(e2[sel2])/np.sqrt(len(e2[sel2]))]
            meane1int += [np.mean(e1int[sel1])]
            meane2int += [np.mean(e2int[sel2])]
            meang1 += [np.mean(g1[sel1])]
            meang2 += [np.mean(g2[sel2])]
            varg1 += [np.std(g1[sel1])/np.sqrt(len(g1[sel1]))]
            varg2 += [np.std(g2[sel2])/np.sqrt(len(g2[sel2]))]
            true_sheared_shape[0] += [np.mean(e1int[sel1]+g1[sel1])]
            true_sheared_shape[1] += [np.mean(e2int[sel2]+g2[sel2])]
            vare1int += [np.std(e1int[sel1])/np.sqrt(len(e1int[sel1]))]
            vare2int += [np.std(e2int[sel2])/np.sqrt(len(e2int[sel2]))]
          #  plt.axvline(bins[i],color='k',linestyle = '--')


        de1 = np.array(meane1int)-np.array(meang1)
        de2 = np.array(meane2int)-np.array(meang2)
        vare1=np.array(vare1)
        vare2=np.array(vare2)
        varg1=np.array(varg1)
        varg2=np.array(varg2)
        vard1 = np.sqrt(vare1*vare1 + varg1*varg1)
        vard2 = np.sqrt(vare2*vare2 + varg2*varg2)
    
        p11, cov11 = np.polyfit(meang1,de1, 1, w=1./vard1/vard1, cov=True, full=False)
        a11,b11 = p11
        p22, cov22 = np.polyfit(meang2,de2, 1, w=1.0/vard2/vard2, cov=True, full=False)
        a22,b22 = p22
    
        print 'Using true sheared shape'
        print 'm11_tss=%f +- %f'%(a11,cov11[0,0]**0.5)
        print 'm22_tss=%f +- %f'%(a22,cov22[0,0]**0.5)
        print 'c11_tss=%f +- %f'%(b11,cov11[1,1]**0.5)
        print 'c22_tss=%f +- %f'%(b22,cov22[1,1]**0.5)
    
        if mode=="show":
            plt.errorbar(true_sheared_shape[0],de1,yerr=vard1, xerr=vare1int, color='m', fmt = 'o', label='$i,j=1,1$')
            plt.errorbar(true_sheared_shape[1],de2,yerr=vard2, xerr=vare2int, color='b', fmt='o', label='$i,j=2,2$')
            plt.plot(true_sheared_shape[0],np.array(meane1int)*(1+a11)+b11, c='m')
            plt.plot(true_sheared_shape[1],np.array(meane2int)*(1+a22)+b22, c='b')
        
            plt.axvline(0.,color='k')
            plt.axhline(0.,color='k')
            plt.legend(loc='upper left')
            plt.xlabel('$<e_i^{int} + g_i>$')
            plt.ylabel('$<e_j>$')
            plt.xlim(bins[0],bins[-1])
        if savedir:
            f = os.path.join(savedir, "%seint-vs-g-truth-plot-auto.png"%bin)
            plt.savefig(f)
            plt.close()
        if mode=='show':
            plt.show()
            plt.close()

    if psf_leakage:
        e1 = res[mask]['e1']
        e2 = res[mask]['e2']

        sel= ((e1psf<1) & (e1psf>-1)) & ((e1<1) & (e1>-1)) & ((e2<1) & (e2>-1)) &  ((e2psf<1) & (e2psf>-1))
        e1psf = e1psf[sel]
        e1 = e1[sel]
        e2 = e2[sel]
        e2psf = e2psf[sel]
        bins = np.linspace(-0.02, 0.02, nbins+1)
        if e1psf.min()>bins[1] or e1psf.max()<bins[-2] or e2psf.min()>bins[1] or e2psf.max()<bins[-2]:
            nbins -=1
        #bins = np.linspace(-0.03, 0.03, nbins)
        if e1psf.min()>bins[1] or e1psf.max()<bins[-2] or e2psf.min()>bins[1] or e2psf.max()<bins[-2]:
            #nbins -=1
            bins = np.linspace(-0.01, 0.01, nbins-1)
        
        meane1=[]; meane2=[]
        vare1=[] ;vare2=[]
        meane1_cross=[]; meane2_cross=[]
        vare1_cross=[] ;vare2_cross=[]
        meane1psf=[] ; meane2psf=[]
        vare1psf=[] ; vare2psf=[]
    
        for i in xrange(len(bins)-1):
            sel1 = (e1psf>bins[i]) & (e1psf<bins[i+1])
            sel2 = (e2psf>bins[i]) & (e2psf<bins[i+1])
            meane1 += [np.mean(e1[sel1])]
            meane2 += [np.mean(e2[sel2])]
            meane1_cross += [np.mean(e1[sel2])]
            meane2_cross += [np.mean(e2[sel1])]
            vare1 += [np.std(e1[sel1])/np.sqrt(len(e1[sel1]))]
            vare2 += [np.std(e2[sel2])/np.sqrt(len(e2[sel2]))]
            vare1_cross += [np.std(e1[sel2])/np.sqrt(len(e1[sel2]))]
            vare2_cross += [np.std(e2[sel1])/np.sqrt(len(e2[sel1]))]
            meane1psf += [np.mean(e1psf[sel1])]
            meane2psf += [np.mean(e2psf[sel2])]
            vare1psf += [np.std(e1psf[sel1])/np.sqrt(len(e1psf[sel1]))]
            vare2psf += [np.std(e2psf[sel2])/np.sqrt(len(e2psf[sel2]))]
            #plt.axvline(bins[i],color='k',linestyle = '--')

        vare1=np.array(vare1)
        vare2 = np.array(vare2)

        p11, cova11 = np.polyfit(e1psf,e1, 1, cov=True, full=False)
        alpha11,c11psf = p11
        p22, cova22 = np.polyfit(e2psf,e2, 1, cov=True, full=False)
        alpha22,c22psf = p22
    
        p12, cova12 = np.polyfit(meane1psf,meane2, 1, cov=True, full=False)
        alpha12,c12psf = p12
        p21, cova21 = np.polyfit(meane2psf,meane1, 1, cov=True, full=False)
        alpha21,c21psf = p21

        print 'PSF correlations:'
        print 'alpha11=%f +- %f'%(alpha11,cova11[0,0]**0.5)
        print 'alpha22=%f +- %f'%(alpha22,cova22[0,0]**0.5)
        print 'alpha12=%f +- %f'%(alpha12,cova12[0,0]**0.5)
        print 'alpha21=%f +- %f'%(alpha21,cova21[0,0]**0.5)
        print 'c11_psf=%f +- %f'%(c11psf,cova11[1,1]**0.5)
        print 'c22_psf=%f +- %f'%(c22psf,cova22[1,1]**0.5)
        print 'c12_psf=%f +- %f'%(c12psf,cova12[1,1]**0.5)
        print 'c21_psf=%f +- %f'%(c21psf,cova21[1,1]**0.5)
    
        if mode=="show":
            plt.errorbar(meane1psf,meane1,yerr=vare1, color='m', fmt = 'o', label='$i,j=1,1$')
            plt.plot(meane1psf, np.array(meane1psf)*alpha11+c11psf, c='m')
            plt.errorbar(meane2psf,meane2,yerr=vare2, color='b', fmt='o', label='$i,j=2,2$')
            plt.plot(meane2psf, np.array(meane2psf)*alpha22+c22psf, c='b')
            plt.axvline(0.,color='k')
            plt.axhline(0.,color='k')
            plt.legend(loc='upper left')
            plt.xlabel('$<e_i^{psf}>$')
            plt.ylabel('$<e_j>$')
            plt.xlim(bins[0],bins[-1])
            plt.ylim(-0.015, 0.015)
        if savedir:
            f = os.path.join(savedir, "%s_epsf-vs-g-auto.png"%bin)
            plt.savefig(f)
        if mode=='show':
            plt.show()
            plt.close()    

        if mode=="show":
            plt.errorbar(meane1psf,meane2_cross,yerr=vare1, color='m', fmt = 'o', label='$i,j=1,2$')
            plt.errorbar(meane2psf,meane1_cross,yerr=vare2, color='b', fmt='o', label='$i,j=2,1$')
            plt.plot(meane1psf, np.array(meane1psf)*alpha12+c12psf, c='m')
            plt.plot(meane2psf, np.array(meane2psf)*alpha21+c21psf, c='b')
        
            plt.axvline(0.,color='k')
            plt.axhline(0.,color='k')
            plt.legend(loc='upper left')
            plt.xlabel('$<e_i^{psf}>$')
            plt.ylabel('$<e_j>$')
            plt.xlim(bins[0],bins[-1])
        if savedir:
            f = os.path.join(savedir, "%s_epsf-vs-g-cross.png"%bin)
            plt.savefig(f)
        if mode=='show':
            plt.show()
            plt.close()
    else:
        alpha11=0.0
        cova11=np.zeros((2,2))
        alpha22=0.0
        cova22=np.zeros((2,2))

        c11psf=0.0
        c22psf=0.

    if return_vals:
        if isinstance(return_vals,bool):
            return (m11, c11, covm11[0,0]**0.5, covm11[1,1]**0.5),(m22, c22, covm22[0,0]**0.5, covm22[1,1]**0.5),(m12,c12,covm12[0,0]**0.5, covm12[1,1]**0.5), (m21,c21,covm21[0,0]**0.5,covm21[1,1]**0.5 ), (alpha11, c11psf, cova11[0,0]**0.5, cova11[1,1]**0.5), (alpha22, c22psf, cova22[0,0]**0.5, cova22[1,1]**0.5), fb
        
        elif return_vals=="alpha":
            return (alpha11+alpha22)/2.0, np.sqrt(cova11[0,0]+cova22[0,0])
        elif return_vals=="alpha1":
            return alpha11, np.sqrt(cova11[0,0])
        elif return_vals=="alpha2":
            return alpha22, np.sqrt(cova22[0,0])

    return 0
    
def histograms_2d(savedir, res):
    """Plot everything measured against everything else measured."""
    os.system("mkdir -p %s/correlations"%savedir)
    for i, t1 in enumerate(res.dtype.names):
        for j, t2 in enumerate(res.dtype.names):
            if (t1 in relevant_parameters) and (t2 in relevant_parameters) and (j>i):
                print "%d %d parameters: %s %s" %(i,j, t1,t2)
                plt.hist2d(res[t1], res[t2],bins=80)
                plt.xlabel(t1)
                plt.ylabel(t2)
                outfile=os.path.join(savedir, "hist2d/hist2d-%s-%s.png"%(t1,t2))
                plt.savefig(outfile)
                plt.close()

def plot_correlations_obs(savedir, res):
    """Plot everything measured against everything else measured."""
    os.system("mkdir -p %s/correlations"%savedir)
    for i, t1 in enumerate(res.dtype.names):
        for j, t2 in enumerate(res.dtype.names):
            if (t1 in relevant_parameters) and (t2 in relevant_parameters) and (j>i):
                print "%d %d parameters: %s %s" %(i,j, t1,t2)
                plt.scatter(res[t1], res[t2])
                plt.xlabel(t1)
                plt.ylabel(t2)
                outfile=os.path.join(savedir, "correlations/%s-%s.png"%(t1,t2))
                plt.savefig(outfile)
                plt.close()

def histograms_1d(savedir, res_sim, res_data=None):
    """Plot histograms of observed quantities for simulation and data"""
    os.system("mkdir -p %s/histograms"%savedir)
    for i, t in enumerate(res_sim.dtype.names):
        if (t in relevant_parameters):
            print "%d %s" %(i, t)
            plt.hist(res_sim[t], bins=80, histtype="step", color="m", normed=1.0, label="simulation")
            if res_data:
                plt.hist(res_sim[t], bins=80, histtype="step", color="m", normed=1.0, label="data")
            plt.xlabel(t)
            outfile=os.path.join(savedir, "histograms/%s.png"%t)
            plt.savefig(outfile)
            plt.close()
  


def extract_weighted_cosmos_sample(catfile, cosmosfile= "/share/des/disc2/image-sims/COSMOS_25.2_training_sample/real_galaxy_catalog_25.2_fits.fits"):
    dat=pyfits.getdata(catfile)
    dist, edges = np.histogram(cat['mag_auto'], bins=np.linspace(16., 26., 50), density=True)
    dat=pyfits.getdata(cosmosfile)
    bin_midpoints = edges[:-1] + np.diff(edges)/2
    cdf = np.cumsum(dist)
    cdf = cdf / cdf[-1]
    values = np.random.rand(5000)
    value_bins = np.searchsorted(cdf, values)
    random_from_cdf = bin_midpoints[value_bins]
    upper, lower=edges[1:], edges[:-1]
    e1=lower[value_bins]
    e2=upper[value_bins]
    lim=16.11, 25.2
    sel=(e1>lim[0]) & (e2<lim[1])
    e1=e1[sel]
    e2=e2[sel]
    ident=[]
    y=[]
    for i in xrange(len(e2)): 
        sel=np.random.randint(0, len(dat['mag_auto'][(dat['mag_auto']>e1[i]) & (dat['mag_auto']<e2[i])]), 1)
        ident+= list(dat['ident'][(dat['mag_auto']>e1[i]) & (dat['mag_auto']<e2[i])][sel])
        y+=list(dat['mag_auto'][(dat['mag_auto']>e1[i]) & (dat['mag_auto']<e2[i])][sel])
        print i
    
    return np.array(ident), np.array(y)    


def quad(a,b):
    return np.sqrt(a*a + b*b)

def gaussian(x,sigma):
    """Return a zero-centred Gaussian of given width"""
    return np.exp(-1. * x * x /2 /sigma /sigma) / np.sqrt(2 * np.pi * sigma * sigma)


def plot_mc_snr(dat=None, comp='', bias='m'):
    if not dat:
        dat=pyfits.getdata('bias_table.fits')

    if isinstance(dat,str):
        dat=pyfits.getdata(dat)
    for i in xrange(len(np.unique(dat['i']))):
        b=dat['%s%s'%(bias,comp)][dat['i']==i]
        eb=dat['err_%s%s'%(bias,comp)][dat['i']==i]
        snr1=dat['vsnr_min'][dat['i']==i]
        snr2=dat['vsnr_max'][dat['i']==i]
        x=(snr1+snr2)/2.
        r1=dat['vrgp_min'][dat['i']==i][0]
        r2=dat['vrgp_max'][dat['i']==i][0]
        plt.errorbar(x,b,yerr=eb, fmt='o', label="$R_{gpp}/R_p=[%2.3f-%2.3f]$"%(r1,r2))
    
    plt.xscale('log')
    plt.ylabel("bias $%s$"%bias)
    if comp:
        plt.ylabel("bias $%s_%d$"%(bias,comp))
    plt.xlabel("SNR")
    plt.legend(loc='lower right')
    if bias=='m':
        plt.ylim(-0.45, 0.15)
    elif bias=='c':
        plt.ylim(-0.004, 0.004)
    plt.axhline(0., color='k')
    plt.show()


def psf_leakage(res, savedir=None, mode="show", equal_number_bins=False):
    mask = np.ones_like(res).astype(bool)

    e1 = res[mask]['e1']
    e2 = res[mask]['e2']
    e1psf = res[mask]['mean_psf_e1_sky']
    e2psf = res[mask]['mean_psf_e2_sky']
    sel= ((e1psf<1) & (e1psf>-1)) & ((e1<1) & (e1>-1)) & ((e2<1) & (e2>-1)) &  ((e2psf<1) & (e2psf>-1))
    e1psf = e1psf[sel]
    e1 = e1[sel]
    e2 = e2[sel]
    e2psf = e2psf[sel]
    nbins = 6

    if not equal_number_bins:
        bins = np.linspace(-0.03, 0.03, nbins+1)
    else:
        bins = find_bin_edges(e1psf, nbins)


    if e1psf.min()>bins[1] or e1psf.max()<bins[-2] or e2psf.min()>bins[1] or e2psf.max()<bins[-2]:
        nbins -=1
    #bins = np.linspace(-0.03, 0.03, nbins)
    if e1psf.min()>bins[1] or e1psf.max()<bins[-2] or e2psf.min()>bins[1] or e2psf.max()<bins[-2]:
        bins = np.linspace(-0.01, 0.01, nbins-1)        
    meane1=[]; meane2=[]
    vare1=[] ;vare2=[]
    meane1_cross=[]; meane2_cross=[]
    vare1_cross=[] ;vare2_cross=[]
    meane1psf=[] ; meane2psf=[]
    vare1psf=[] ; vare2psf=[]    
    for i in xrange(len(bins)-1):
        sel1 = (e1psf>bins[i]) & (e1psf<bins[i+1])
        sel2 = (e2psf>bins[i]) & (e2psf<bins[i+1])
        meane1 += [np.mean(e1[sel1])]
        meane2 += [np.mean(e2[sel2])]
        meane1_cross += [np.mean(e1[sel2])]
        meane2_cross += [np.mean(e2[sel1])]
        vare1 += [np.std(e1[sel1])/np.sqrt(len(e1[sel1]))]
        vare2 += [np.std(e2[sel2])/np.sqrt(len(e2[sel2]))]
        vare1_cross += [np.std(e1[sel2])/np.sqrt(len(e1[sel2]))]
        vare2_cross += [np.std(e2[sel1])/np.sqrt(len(e2[sel1]))]
        meane1psf += [np.mean(e1psf[sel1])]
        meane2psf += [np.mean(e2psf[sel2])]
        vare1psf += [np.std(e1psf[sel1])/np.sqrt(len(e1psf[sel1]))]
        vare2psf += [np.std(e2psf[sel2])/np.sqrt(len(e2psf[sel2]))]
        plt.axvline(bins[i],color='k',linestyle = '--')

    vare1=np.array(vare1)
    vare2 = np.array(vare2)
    p11, cova11 = np.polyfit(meane1psf,meane1, 1, cov=True, w=1./vare1/vare1,full=False)
    alpha11,c11psf = p11
    p22, cova22 = np.polyfit(meane2psf,meane2, 1, cov=True, w=1./vare2/vare2,full=False)
    alpha22,c22psf = p22    
    p12, cova12 = np.polyfit(meane1psf,meane2, 1, cov=True, w=1./vare1/vare1,full=False)
    alpha12,c12psf = p12
    p21, cova21 = np.polyfit(meane2psf,meane1, 1, cov=True, w=1./vare2/vare2,full=False)
    alpha21,c21psf = p21

    print 'PSF correlations:'
    print 'alpha11=%f +- %f'%(alpha11,cova11[0,0]**0.5)
    print 'alpha22=%f +- %f'%(alpha22,cova22[0,0]**0.5)
    print 'alpha12=%f +- %f'%(alpha12,cova12[0,0]**0.5)
    print 'alpha21=%f +- %f'%(alpha21,cova21[0,0]**0.5)
    print 'c11_psf=%f +- %f'%(c11psf,cova11[1,1]**0.5)
    print 'c22_psf=%f +- %f'%(c22psf,cova22[1,1]**0.5)
    print 'c12_psf=%f +- %f'%(c12psf,cova12[1,1]**0.5)
    print 'c21_psf=%f +- %f'%(c21psf,cova21[1,1]**0.5)    
    plt.errorbar(meane1psf,meane1,yerr=vare1, color='m', fmt = 'o', label='$i,j=1,1$')
    plt.plot(meane1psf, np.array(meane1psf)*alpha11+c11psf, c='m')
    plt.errorbar(meane2psf,meane2,yerr=vare2, color='b', fmt='o', label='$i,j=2,2$')
    plt.plot(meane2psf, np.array(meane2psf)*alpha22+c22psf, c='b')
    plt.axvline(0.,color='k')
    plt.axhline(0.,color='k')
    plt.legend(loc='upper left')
    plt.xlabel('$<e_i^{psf}>$')
    plt.ylabel('$<e_j>$')
    plt.xlim(-0.02,0.02)
    plt.ylim(-0.005, 0.005)

    if savedir:
        f = os.path.join(savedir, "epsf-vs-g-auto.png")
        plt.savefig(f)

    if mode=='show':
        plt.show()
        plt.close()    
    plt.errorbar(meane1psf,meane2_cross,yerr=vare1, color='m', fmt = 'o', label='$i,j=1,2$')
    plt.errorbar(meane2psf,meane1_cross,yerr=vare2, color='b', fmt='o', label='$i,j=2,1$')
    plt.plot(meane1psf, np.array(meane1psf)*alpha12+c12psf, c='m')
    plt.plot(meane2psf, np.array(meane2psf)*alpha21+c21psf, c='b')    
    plt.axvline(0.,color='k')
    plt.axhline(0.,color='k')
    plt.legend(loc='upper left')
    plt.xlabel('$<e_i^{psf}>$')
    plt.ylabel('$<e_j>$')
    plt.xlim(-0.02,0.02)
    plt.ylim(-0.005, 0.005)

    if savedir:
        f = os.path.join(savedir, "%s_epsf-vs-g-cross.png"%bin)
        plt.savefig(f)
    if mode=='show':
        plt.show()
        plt.close()

    return (alpha11, cova11[0,0]), (alpha22, cova22[0,0]), (c11psf, cova11[1,1]), (c22psf, cova22[1,1]) 

def bootstrap_error(nsubsamples, full_cat, operation, additional_args=None, additional_argvals=None, method="split", sample_frac=0.1, columns_needed=[]):
    """Generic function which takes an array and splits into nsubsamples parts. Then apply an operation
       to each, and find the rms deviation about the mean in the resulting quantity."""

    methods=["sky_coord", "split", "random", "cosmos"]
    if method not in methods:
        raise ValueError("Unrecognised method: %s"%method)
    else:
        print "Generating bootstrap errorbars: %s"%method
        print "Using %d samples"%nsubsamples
    resampled = []
    if isinstance(full_cat, tuple) and (len(full_cat)==2):
        if len(columns_needed)>0:
            print "Using selected columns"
            dtx= [ ( n, str(full_cat[0].dtype[n]) ) for n in full_cat[0].dtype.names if n in columns_needed] 
            dty= [ ( n, str(full_cat[1].dtype[n]) ) for n in full_cat[1].dtype.names if n in columns_needed]
            dt = np.dtype(dtx + dty)
            data = np.zeros(full_cat[0].size, dtype=dt)
            xdata = full_cat[0]
            ydata = full_cat[1]
            
            for col in np.dtype(dtx).names: 
                data[col] = full_cat[0][col]
            for col in np.dtype(dty).names:
                data[col] = full_cat[1][col]
            merged_table = True
        else:
            xdata = full_cat[0]
            ydata = full_cat[1]
            merged_table = False
        separate_xy = True
       
    else:
        separate_xy = False
        ydata = full_cat
    if method=="split":
        resample_length=len(ydata)/nsubsamples
    else:
        ntot = ydata.size
        resample_length = int(sample_frac * ntot)

    if (not method=="sky_coord"):
        bootstrap_edges=(np.arange(nsubsamples+1)*resample_length).astype(int)
    else:
        # Split the footprint into a set of ra rows
        ndir = int(np.sqrt(nsubsamples))
        ra_bins = find_bin_edges(ydata["ra"], ndir)
        ra_bins = zip(ra_bins[:-1], ra_bins[1:] )
        print "Splitting data into %dx%d bootstrap patches"%(ndir,ndir)

        ira = 0
        idec= 0

    for i in xrange(nsubsamples):
        print "    %d"%(i+1)

        if (method=="split"):
            b_low = bootstrap_edges[i]
            b_high = bootstrap_edges[i+1]
            indices = np.arange(b_low, b_high, 1)
        elif (method=="sky_coord"):

            if ira==ndir: continue

            # Galaxies in this ra bin
            sel_ra = (ydata["ra"]>ra_bins[ira][0]) & (ydata["ra"]<ra_bins[ira][1])

            dec_bins =  find_bin_edges(ydata["dec"][sel_ra], ndir)
            dec_bins = zip(dec_bins[:-1], dec_bins[1:] )
            # Galaxies in this dec bin
            sel_dec = (ydata["dec"]>dec_bins[idec][0]) & (ydata["dec"]<dec_bins[idec][1])

            # Cut this region out and repeat the calculation using only the galaxies outside it
            indices = np.invert(sel_ra & sel_dec)

            if idec==ndir-1:
                idec=0
                ira+=1
            else:
                idec+=1

        elif (method=="random"):
            indices = np.random.rand(ntot)<sample_frac #np.random.choice(ntot-1,resample_length)
        elif (method=="cosmos"):
            print "Subdividing catalogue by COSMOS identifier"
            cosmos_ids = np.unique(xdata["cosmos_ident"])
            selected_ids = np.random.choice(cosmos_ids, cosmos_ids.size/2, replace=False)
            indices = np.in1d(xdata["cosmos_ident"], selected_ids)
        #print "bootstrap subsample %d (%d-%d)"%(i+1, b_high, b_low)
        if additional_args is None:
            if not separate_xy:
                derived_quantity = operation( full_cat[indices])
            else:
                derived_quantity = operation( xdata[indices], ydata[indices])
        else:
            kwargs = {}
            for name,val in zip(additional_args, additional_argvals):
                kwargs[name] = val



            if not separate_xy:
                derived_quantity = operation(full_cat[indices], **kwargs)
            else:
                if not merged_table:
                    derived_quantity = operation(xdata[indices], ydata[indices], **kwargs)
                else:
                    dat = data[indices]
                    derived_quantity = operation(dat, dat, **kwargs)

            #print derived_quantity


        resampled.append(derived_quantity)

    fmean = np.mean(resampled)
    dev = (np.array(resampled)-fmean) * (np.array(resampled)-fmean)

    return np.sqrt(np.mean(dev))

def find_bin_edges(x,nbins,w=None):
    """
    For an array x, returns the boundaries of nbins equal (possibly weighted by w) bins.
    """

    if w is None:
      xs=np.sort(x)
      r=np.linspace(0.,1.,nbins+1.)*(len(x)-1)
      return xs[r.astype(int)]

    fail=False
    ww=np.sum(w)/nbins
    i=np.argsort(x)
    k=np.linspace(0.,1.,nbins+1.)*(len(x)-1)
    k=k.astype(int)
    r=np.zeros((nbins+1))
    ist=0
    for j in xrange(1,nbins):
      # print k[j],r[j-1]
      if k[j]<r[j-1]:
        print 'Random weight approx. failed - attempting brute force approach'
        fail=True
        break
      w0=np.sum(w[i[ist:k[j]]])
      if w0<=ww:
        for l in xrange(k[j],len(x)):
          w0+=w[i[l]]
          if w0>ww:
            r[j]=x[i[l]]
            ist=l
            break
      else:
        for l in xrange(k[j],0,-1):
          w0-=w[i[l]]
          if w0<ww:
            r[j]=x[i[l]]
            ist=l
            break

    if fail:

      ist=np.zeros((nbins+1))
      ist[0]=0
      for j in xrange(1,nbins):
        wsum=0.
        for k in xrange(ist[j-1].astype(int),len(x)):
          wsum+=w[i[k]]
          if wsum>ww:
            r[j]=x[i[k-1]]
            ist[j]=k
            break

    r[0]=x[i[0]]
    r[-1]=x[i[-1]]

    return r


def im3shape_weights_grid(data, bins_from_table=True, table=None, filename="/home/samuroff/shear_pipeline/im3shape-weights_column.fits", sbins=15, rbins=15, binning="equal_number", simdat=None):

    if len(data.dtype.names)>3:
        print "Extracting columns required"
        placeholder=data
        data = np.zeros(data.size, dtype=[("snr", float ), ("mean_rgpp_rp", float), ("e1", float)])
        for col in data.dtype.names:
            data[col] = placeholder[col]

    
    if bins_from_table:
        bt = fi.FITS(table)[1].read()
        rgpp = (np.unique(bt["rgp_lower"])+np.unique(bt["rgp_upper"]))/2
        nbins = rgpp.size

        print "Will evaluate weights on grid from %s"%table
    elif (binning is "equal_number"):
        bt=np.zeros(rbins, dtype=[("rgp_lower", float), ("rgp_upper", float)])
        bins = find_bin_edges(data["mean_rgpp_rp"], rbins)
        bt["rgp_lower"] = bins[:-1]
        bt["rgp_upper"] = bins[1:]
    else:
        print "Using log bins in rgpp"
        bt=np.zeros(rbins, dtype=[("rgp_lower", float), ("rgp_upper", float)])
        bt["rgp_lower"] = np.logspace(np.log10(1.13),np.log10(3),rbins+1)[:-1]
        bt["rgp_upper"] = np.logspace(np.log10(1.13),np.log10(3),rbins+1)[1:]
    
    wt_vec=[]
    for i, (rgpp_lower, rgpp_upper) in enumerate(zip(np.unique(bt["rgp_lower"]), np.unique(bt["rgp_upper"]))):

        row_data = data[ (data["mean_rgpp_rp"]>rgpp_lower) & (data["mean_rgpp_rp"]<rgpp_upper) ]

        if simdat is not None:
            sim_row_data = simdat[ (simdat["mean_rgpp_rp"]>rgpp_lower) & (simdat["mean_rgpp_rp"]<rgpp_upper) ]

        if bins_from_table:
            select_row = (bt["i"]==i)
            snr_edges0 = bt["snr_lower"][select_row]
            snr_edges1 = bt["snr_upper"][select_row]
        elif (binning is "equal-number"):
            snr_edges = find_bin_edges(row_data["snr"], sbins)
            snr_edges0 = snr_edges[:-1]
            snr_edges1 = snr_edges[1:]
        else:
            snr_edges0 = np.logspace(np.log10(12),np.log10(200),sbins+1)[:-1]
            snr_edges1 = np.logspace(np.log10(12),np.log10(200),sbins+1)[1:]

        for j, (snr_lower, snr_upper) in enumerate(zip(snr_edges0, snr_edges1)):

            cell_sample = (row_data["snr"]>snr_lower) &  (row_data["snr"]<snr_upper)

            wt = compute_im3shape_weight(row_data["e1"][cell_sample], verbose=False)
            ngal  = row_data["e1"][cell_sample].size

            if simdat is not None:
                sim_cell_sample= (sim_row_data["snr"]>snr_lower) &  (sim_row_data["snr"]<snr_upper)
                wt_sim = compute_im3shape_weight(sim_row_data["e1"][sim_cell_sample], verbose=False)
                nsim = sim_row_data["e1"][sim_cell_sample].size
            else:
                wt_sim = 0
                nsim = 0

            print "%d,%d [%3.2f,%3.2f], [%3.2f,%3.2f] %d"%(i,j,rgpp_lower, rgpp_upper, snr_lower, snr_upper, ngal)

            wt_vec.append([i, j, rgpp_lower, rgpp_upper, snr_lower, snr_upper, ngal, wt, nsim, wt_sim ])

    if filename is not None:
        out_fits = fi.FITS(filename, "rw")

        out_table = np.zeros(len(wt_vec), dtype=[("i_rgpp", int),("i_snr", int), ("rgpp_lower", float),("rgpp_upper", float), ("snr_lower", float),("snr_upper", float), ("ngal", float), ("inverse_weight", float), ("nsim", float), ("simulation_inverse_weight", float)])
        for i, col in enumerate(out_table.dtype.names):
            print "Writing column %s"%col
            out_table[col] = np.array(wt_vec).T[i]

        out_fits.write(out_table, clobber=True)
        out_fits.close()
    return out_table




def compute_im3shape_weight(e1, verbose=True):
    """Return either the direct measurement of the variance, or the width of the best fit Gaussian distribution"""
    sigma_direct = e1.std()

    data, bins = np.histogram(e1, normed=1, bins=50)
    x = (bins[1:]+bins[:-1])/2
    param, cov = optimise.curve_fit(gaussian, x, data, [0.25])
    sigma_fit = param[0]

    if (not np.isfinite(sigma_fit)) and (np.isfinite(sigma_direct)):
        if verbose:
            print "WARNING fit variance is not finite: %2.4f"%simga_direct
        return sigma_direct

    elif (np.isfinite(sigma_fit)) and (not np.isfinite(sigma_direct)):
        if verbose:
            print "WARNING direct variance estimate is not finite: %2.4f"%sigma_fit
        return sigma_fit
    elif (np.isfinite(sigma_fit)) and (np.isfinite(sigma_direct)):
        if verbose:
            print "Using the maximum of direct : %2.4f, best fit : %2.4f"%(sigma_direct,sigma_fit)
        return max(sigma_direct,sigma_fit)

    else:
        print "WARNING: neither variance estimate is finite. Using default (0.25)"
        return 0.25

def interpolate_weights_grid(weights_grid, target_data, smoothing=1.0, outdir=".", outfile="im3shape_weights.fits"):
    from scipy.interpolate import Rbf
    # Set up the interpolator with some fixed scale factor
    print "Setting up variance interpolator"
    rgpp = np.sqrt(weights_grid["rgpp_lower"] * weights_grid["rgpp_upper"])
    snr = np.sqrt(weights_grid["snr_lower"] * weights_grid["snr_upper"])
    fx = np.log10(snr).max()
    fy = np.log10(rgpp).max()

    interpolator = Rbf(np.log10(snr)/fx, np.log10(rgpp)/fy, weights_grid["inverse_weight"], smooth=smoothing, function="multiquadric")

    import pdb ; pdb.set_trace()

    print "Performing interpolation"
    ngal = target_data["snr"].size
    if ngal<10e6:
        interpolated_sigma = interpolator( np.log10(target_data["snr"])/fx, np.log10(target_data["mean_rgpp_rp"])/fy )
    else:
        print "Splitting array"
        ngal = target_data["snr"].size
        nchunk = ngal/4
        interpolated_sigma=[]
        for i in xrange(4):
            print i
            interpolated_sigma.append( interpolator( np.log10(target_data["snr"][i*nchunk:(i+1)*nchunk])/fx, np.log10(target_data["mean_rgpp_rp"][i*nchunk:(i+1)*nchunk])/fy ) )
        if ngal!=4*nchunk:
            interpolated_sigma.append( interpolator( np.log10(target_data["snr"][(i+1)*nchunk:])/fx, np.log10(target_data["mean_rgpp_rp"][(i+1)*nchunk:])/fy ) )

        interpolated_sigma = np.concatenate(interpolated_sigma)

    print "Setting up output arrays"
    out = np.empty(target_data.size, dtype=[("coadd_objects_id", int), ("weight", float)])
    out["coadd_objects_id"] = target_data["coadd_objects_id"] 
    out["weight"] = 1./(interpolated_sigma * interpolated_sigma)

    print "Saving output to %s/%s"%(outdir, outfile)
    out_fits = fi.FITS("%s/%s"%(outdir, outfile), "rw")
    out_fits.write(out)
    out_fits[-1].write_key("EXTNAME", "i3s_weights_col")
    out_fits.close()
    print "Done"











def get_bias(xdata, catalogue, external_calibration_col=None, apply_calibration=False, nbins=5, weights=None, xlim=(-1,1), ellipticity_name="e", names=["m","c"], visual=False, binning="equal_number", silent=False, mask=None):

    g1 = xdata['true_g1']
    g2 = xdata['true_g2']

    
    sel = (g1>xlim[0]) & (g1<xlim[1]) & (g2>xlim[0]) & (g2<xlim[1])
    if apply_calibration and (external_calibration_col is None):
        try:
            mask_nans = np.isfinite(catalogue["m"]) & np.isfinite(catalogue["c1"]) & np.isfinite(catalogue["c2"]) 
        except:
            mask_nans = np.isfinite(catalogue["m"])
        sel = sel & mask_nans
    elif (external_calibration_col is not None):
        mask_nans = np.isfinite(external_calibration_col["m"]) & np.isfinite(external_calibration_col["c1"]) & np.isfinite(external_calibration_col["c2"]) 
        sel = sel & mask_nans
    if mask is not None:
        sel = sel & mask
    g1 = g1[sel]
    g2 = g2[sel]

    try:
        e1 = catalogue["%s1"%ellipticity_name][sel]
        e2 = catalogue["%s2"%ellipticity_name][sel]
    except:
        if ellipticity_name is "true_sheared_e":
            e1 = catalogue["intrinsic_e1"][sel] + g1
            e2 = catalogue["intrinsic_e2"][sel] + g2

    if weights is not None:
        print "Will use weights provided"
        w = weights[sel]
    else:
        w = np.ones_like(e1)

    if apply_calibration:
        print "Applying calibration columns"
        if external_calibration_col is not None:
            print "Using external calibration column"
            c1 = external_calibration_col["c1"][sel]
            c2 = external_calibration_col["c2"][sel]
            m = external_calibration_col["m"][sel]
        elif "c1" in catalogue.dtype.names:
            m = catalogue["m"][sel]
            c1 = catalogue["c1"][sel]
            c2 = catalogue["c2"][sel]
        else:
            print "No additive column found"
            c1 = np.zeros_like(e1)
            c2 = np.zeros_like(e1)
    else:
        m = np.zeros_like(e1)
        c1 = np.zeros_like(e1)
        c2 = np.zeros_like(e1)


    if isinstance(binning,str):
        if binning=="equal_number":
            bins = find_bin_edges(g1, nbins)
    else:
        bins = binning

    x = (bins[1:]+bins[:-1])/2.0

    y11,y22,y12,y21=[],[],[],[]
    variance_y11,variance_y22,variance_y12,variance_y21=[],[],[],[]
    
    for i, bounds in enumerate(zip(bins[:-1],bins[1:])):
        lower = bounds[0]
        upper = bounds[1]
        sel1 = (g1>lower) & (g1<upper)
        sel2 = (g2>lower) & (g2<upper)

        y11.append(np.sum(w[sel1]*(e1[sel1]-c1[sel1])) / np.sum(w[sel1]*(1+m[sel1])))
        variance_y11.append( compute_im3shape_weight(e1[sel1]-g1[sel1], verbose=False) / (e1[sel1].size**0.5) )

        y22.append(np.sum(w[sel2]*(e2[sel2]-c2[sel2])) / np.sum(w[sel2]*(1+m[sel2])))
        variance_y22.append( compute_im3shape_weight(e2[sel2]-g2[sel2], verbose=False) / (e2[sel2].size**0.5) )

        y12.append(np.sum(w[sel2]*(e1[sel2]-c1[sel2])) / np.sum(w[sel2]*(1+m[sel2])))
        variance_y12.append( compute_im3shape_weight(e2[sel1]-g2[sel1], verbose=False) / (e2[sel1].size**0.5) )
       
        y21.append(np.sum(w[sel1]*(e2[sel1]-c2[sel1])) / np.sum(w[sel1]*(1+m[sel1])))
        variance_y21.append( compute_im3shape_weight(e1[sel2]-g1[sel2], verbose=False) / (e1[sel2].size**0.5) )

    d1 = np.array(y11)-x
    d2 = np.array(y22)-x

    # Four linear fits to do here

    # m11, c11
    try:
        p11, cov11 = optimise.curve_fit(fline, x, d1, sigma=np.array(variance_y11), p0=[0.05,1e-7])
    except: import pdb ; pdb.set_trace()
    if not np.isfinite(cov11).all():
        p11,cov11 = stabilised_fit(x, d1, w11)
    m11,c11 = p11  

    # m22, c22
    p22, cov22 = optimise.curve_fit(fline, x, d2, sigma=np.array(variance_y22), p0=[0.05,1e-7])
    if not np.isfinite(cov22).all():
        p22,cov22 = stabilised_fit(x,d2,w22)
    m22,c22 = p22

    # m12, c12
    p12, cov12 = optimise.curve_fit(fline, x, y12, sigma=np.array(variance_y12), p0=[0.05,1e-7])
    if not np.isfinite(cov12).all():
        p12,cov12 = stabilised_fit(x,y12,1.0/np.array(variance_y12)/np.array(variance_y12))
    m12,c12 = p12

    # m21, c21
    p21, cov21 = optimise.curve_fit(fline, x, y21, sigma=np.array(variance_y21), p0=[0.05,1e-7])
    if not np.isfinite(cov21).all():
        p21,cov21 = stabilised_fit(x,y21,1.0/np.array(variance_y21)/np.array(variance_y21))
    m21,c21 = p21

    if visual:
        import pylab as plt
        plt.subplots_adjust(left=0.15,right=0.95)
        plt.errorbar(x, d1, yerr=np.array(variance_y11), fmt="o", color="purple", label="$e_1$")
        plt.errorbar(x, d2, yerr=np.array(variance_y22), fmt="D", color="steelblue", label="$e_2$")
        plt.plot(x, x*m11+c11, "-", lw=2.5, color="purple")
        plt.plot(x, x*m22+c22, "-", lw=2.5, color="steelblue")
        plt.xlabel("Input Shear $g_i$")
        plt.ylabel("Residual Shear $<e_i> - <g_i>$")
        plt.legend(loc="lower left")

    if not silent:
        print 'm11=%f +- %f'%(m11,cov11[0,0]**0.5)
        print 'm22=%f +- %f'%(m22,cov22[0,0]**0.5)
        print 'm12=%f +- %f'%(m12,cov12[0,0]**0.5)
        print 'm21=%f +- %f'%(m21,cov21[0,0]**0.5)
        print 'c11=%f +- %f'%(c11,cov11[1,1]**0.5)
        print 'c22=%f +- %f'%(c22,cov22[1,1]**0.5)
        print 'c12=%f +- %f'%(c12,cov12[1,1]**0.5)
        print 'c21=%f +- %f'%(c21,cov21[1,1]**0.5)

    m = (m11+m22)/2
    c = (c11+c22)/2
    error_m = np.sqrt(cov11[0,0] + cov22[0,0])
    error_c = np.sqrt(cov11[1,1] + cov22[1,1])

    biases_dict = {"m":(m,error_m), "c":(c,error_c), "m11":(m11,cov11[0,0]**0.5), "c11":(c11,cov11[1,1]**0.5), "m22":(m22,cov22[0,0]**0.5), "c22":(c22,cov22[1,1]**0.5), "m12":(m12,cov12[0,0]**0.5), "c12":(c12,cov12[1,1]**0.5), "m21":(m21,cov21[0,0]**0.5), "c21":(c21,cov21[1,1]**0.5)}

    # Decide how to package the results
    # Either as a dictionary, including the fit covariances
    # Or as a single value (probably for when using bootstrap errors)
    out={}
    if not isinstance(names, str):
        for name in names:
            out[name] = biases_dict[name]
    else:
        out = biases_dict[names][0]

    return out



def setup_resampling_function(source, nbin):
    """Set up a function which takes a set of bin coefficients and
       resamples a predefined histogram with predefined bin edges."""

    def resampling_function(*args):
        """ Apply a set of bin coefficients to resample a histogram
            The only parameters this function should take is a list
            of cofficients and the bin edges.
            The idea here is that once defined, we can fit this 
            function to an an arbitrary target distribution."""
        bins = args[0]
        bin_fractions = args[1:nbin+1] 
        pq = np.histogram(source, bins=bins)[0]
        return (pq * np.array(bin_fractions) ).astype(int)

    return resampling_function


def get_selection_to_match(target, unweighted, nbins=40, xlim=(None,None), existing_weights=None):
    """Return a selection mask which will resample one distribution to look like another."""

    if xlim[0] is None:
        dx=(unweighted.max()-unweighted.min())/nbins
        xlim=( unweighted.min()-dx, unweighted.max()+dx)

    if target.size>unweighted.size:
        target = target[:unweighted.size]

    n_uw, bins_uw = np.histogram(unweighted, bins=nbins, range=xlim, weights=existing_weights)
    n_tar, bins_tar = np.histogram(target, bins=nbins, range=xlim)

    # Set up a blank selection function to fill
    mask = np.ones_like(unweighted).astype(bool)

    centres=(bins_uw[1:]+bins_uw[:-1])/2

    frac = n_tar*1.0/n_uw / (n_tar*1.0/n_uw).max() 

    ntot=unweighted.size

    for i, edges in enumerate(zip(bins_uw[:-1], bins_uw[1:])):
        lower = edges[0]
        upper = edges[1]
        print "%d [%2.3f-%2.3f]"%(i,lower,upper)

        select_subset = (unweighted<upper) & (unweighted>lower)
        nbin = mask[select_subset].size
        mask_subset = np.ones(nbin)

        bin_frac = frac[i]

        nremove = nbin - bin_frac*nbin
        nremove = int(nremove)

        print "Selecting %d/%d (%2.3f pc)"%(nbin-nremove,nbin, 100.*(1-nremove*1.0/nbin))

        indices_to_keep = np.random.choice(range(nbin), nremove, replace=False)
        mask_subset[indices_to_keep] = 0

        mask[select_subset] = mask_subset.astype(bool)
   

    return mask


def get_selection_to_match_3d(target, unweighted, nbins=40, xlim=(None,None), existing_weights=None):
    """Return a selection mask which will resample one distribution to look like another."""



    if target[0].size>unweighted[0].size:
        target = [t[:unweighted[0].size] for t in target]

    n_uw, bins_uw = np.histogramdd(unweighted, bins=nbins)
    n_tar, bins_tar = np.histogramdd(target, bins=nbins)

    # Set up a blank selection function to fill
    mask = np.ones_like(unweighted[0]).astype(bool)

    frac = n_tar*1.0/n_uw
    frac[np.invert(np.isfinite(frac))]=0
    frac/= frac.max() 

    ndim = len(bins_tar)

    for i, edges0 in enumerate(zip(bins_uw[0][:-1], bins_uw[0][1:])):
        for j, edges1 in enumerate(zip(bins_uw[1][:-1], bins_uw[1][1:])):
            for k, edges2 in enumerate(zip(bins_uw[2][:-1], bins_uw[2][1:])):

                print "%d %d %d [%2.3f-%2.3f] [%2.3f-%2.3f] [%2.3f-%2.3f]"%(i,j,k, edges0[0], edges0[1], edges1[0], edges1[1], edges2[0], edges2[1])

                select_subset = (unweighted[0]<edges0[1]) & (unweighted[0]>edges0[0]) & (unweighted[1]<edges1[1]) & (unweighted[1]>edges1[0]) & (unweighted[2]<edges2[1]) & (unweighted[2]>edges2[0])
                nbin = mask[select_subset].size
                if nbin==0: continue
                mask_subset = np.ones(nbin)

                bin_frac = frac[i,j,k]

                nremove = nbin - bin_frac*nbin
                nremove = int(nremove)

                print "Selecting %d/%d (%2.3f pc)"%(nbin-nremove,nbin, 100.*(1-nremove*1.0/nbin))

                indices_to_keep = np.random.choice(range(nbin), nremove, replace=False)
                mask_subset[indices_to_keep] = 0

                mask[select_subset] = mask_subset.astype(bool)
   

    return mask

def get_weights_to_match(target, unweighted, nbins=60, xlim=(None,None), existing_weights=None):
    """Return a set of weights which will force one distribution to look like another."""

    if xlim[0] is None:
        dx=(unweighted.max()-unweighted.min())/nbins
        xlim=( unweighted.min()-dx, unweighted.max()+dx)

    n_uw, bins_uw = np.histogram(unweighted, bins=nbins, normed=1, range=xlim, weights=existing_weights)
    n_tar, bins_tar = np.histogram(target, bins=nbins, normed=1, range=xlim)

    import scipy.interpolate
    x_tar = (bins_tar[1:]+bins_tar[:-1])/2.
    p_tar = scipy.interpolate.interp1d(x_tar, n_tar)

    x_uw = (bins_uw[1:]+bins_uw[:-1])/2.
    p_uw = scipy.interpolate.interp1d(x_uw, n_uw)

    return p_tar(unweighted)/p_uw(unweighted)

def get_2d_wts(target, source, nbins=60, xlim=(None,None), ylim=(None,None), verbose=False):
    print "Constructing histograms."
    if (xlim[0] is not None) or (ylim[0] is not None):
        rge=(xlim,ylim)
    else:
        rge=None
    n_tar, x_tar, y_tar = np.histogram2d(target[0], target[1], bins=nbins, normed=False, range=rge)
    wts = np.zeros(source[0].size)


    if (target[1].size!=source[1].size):
        trunc = min(target[1].size,source[1].size)
        if verbose: 
            print "Enforcing equal sized arrays. Size %d"%trunc

        target_pool = [[],[]]
        source_pool = [[],[]]     
        target_pool[0] = target[0][:trunc]
        target_pool[1] = target[1][:trunc]
        source_pool[0] = source[0][:trunc]
        source_pool[1] = source[1][:trunc]
    else:
        target_pool = target
        source_pool = source

    for i, edges in enumerate(zip(x_tar[:-1], x_tar[1:])):
        sel_t = (target_pool[0]>edges[0]) & (target_pool[0]<edges[1])
        sel_s = (source_pool[0]>edges[0]) & (source_pool[0]<edges[1])
        tg = target_pool[1][sel_t]
        sc = source_pool[1][sel_s]

        for j, edges_inner in enumerate(zip(y_tar[:-1], y_tar[1:])):
            if verbose:
                print "Processing node (%d %d) edges : [%3.3f - %3.3f] [%3.3f - %3.3f]"%(i+1, j+1, edges[0], edges[1], edges_inner[0], edges_inner[1])

            target_inner = (tg>edges_inner[0]) & (tg<edges_inner[1])
            source_inner = (sc>edges_inner[0]) & (sc<edges_inner[1])
           
            try:
                select = (source[1]>edges_inner[0]) & (source[1]<edges_inner[1]) & (source[0]>edges[0]) & (source[0]<edges[1])
                wts[select] = tg[target_inner].size * 1.0 / sc[source_inner].size 
            except:
                select = (source[1]>edges_inner[0]) & (source[1]<edges_inner[1]) & (source[0]>edges[0]) & (source[0]<edges[1])
                wts[select] = np.nan

    wts[np.invert(np.isfinite(wts))]=0

    return wts



def get_nd_wts(target, source, nbins=60, verbose=False):
    print "Constructing histograms."

    # work out the number of dimensions required
    # Define the binning in each
    bounds=[]
    bins = []
    edges=[]
    for dimension in source:
        bounds.append((dimension.min(), dimension.max()))
        bins.append(np.linspace(0, nbins-1, nbins).astype(int))
        bin_centres=np.linspace(bounds[-1][0], bounds[-1][1], nbins+1)
        edges.append(zip( bin_centres[:-1], bin_centres[1:]))
    ndim = len(bounds)
    if verbose:
        print "Found %d dimensions"%ndim
    bins_all = np.meshgrid(*bins)

    wts = np.zeros(source[0].size)

    # Truncate the arrays to ensure the comparison is meaningful
    if (target[1].size!=source[1].size):
        trunc = min(target[1].size,source[1].size)
        if verbose: 
            print "Enforcing equal sized arrays. Size %d"%trunc
        target_pool = []
        source_pool = []     
        for i in xrange(ndim):
            target_pool.append( target[i][:trunc])
            source_pool.append(source[i][:trunc])     
    else:
        target_pool = target
        source_pool = source

    # Now loop over the ndim dimensional grid
    ntot = nbins**ndim

    for i in xrange(ntot):
        if verbose:
                print "Processing node %d/%d edges : "%(i+1, ntot),
        data_source = source_pool
        data_target = target_pool
        selection = np.ones_like(wts).astype(bool)
        for j, dimension in enumerate(zip(bins_all, edges)):
            bin_index = dimension[0].flatten()[i]
            lower,upper = dimension[1][bin_index][0], dimension[1][bin_index][1]

            if verbose: print "[%3.3f - %3.3f]"%(lower,upper),
            selection = selection & (source[j]<upper) & (source[j]>lower)
            selection_source = (data_source[j]<upper) & (data_source[j]>lower)
            data_source = [ dat[selection_source] for dat in data_source]
            selection_target = (data_target[j]<upper) & (data_target[j]>lower)
            data_target = [ dat[selection_target] for dat in data_target]

        try:
            wts[selection] = data_target[0].size * 1.0 / data_source[0].size 
            if verbose: print ""

        except:
            wts[selection] = 0.0
            if verbose: print 0


    wts[np.invert(np.isfinite(wts))]=0

    np.savetxt("weights_column.txt", wts)

    return wts


def get_weights_surface(target, unweighted, nbins=90, xlim=(None,None), ylim=(None,None), split=False, do_diagnostic_plot=False):
    """Return a 2d weight grid which will force one distribution to look like another."""

    print "Constructing histograms."
    n_uw, x_uw , y_uw= np.histogram2d(unweighted[0], unweighted[1], bins=nbins, normed=False, range=(xlim,ylim))
    n_tar, x_tar, y_tar = np.histogram2d(target[0], target[1], bins=nbins, normed=False, range=(xlim,ylim))


    if xlim[1] is not None:
        fx = np.log10(xlim[1])
        fy = np.log10(ylim[1])
    else:
        fx = np.log10(x_tar).max()
        fy = np.log10(y_tar).max()

    import scipy.interpolate
    print "Setting up interpolators...",
    print "1",
    x_tar = (x_tar[1:]+x_tar[:-1])/2.
    y_tar = (y_tar[1:]+y_tar[:-1])/2.
    xx,yy=np.meshgrid(x_tar,y_tar)
    norm = np.trapz([np.trapz(n, y_tar) for n in n_tar], x_tar)
    p_tar = setup_rbf_interpolator(xx.flatten(), yy.flatten(), n_tar.flatten(), logx=True, logy=True, scale=[fx,fy])

    
    print "2"
    x_uw = (x_uw[1:]+x_uw[:-1])/2.
    y_uw = (y_uw[1:]+y_uw[:-1])/2.
    xx0,yy0 = np.meshgrid(x_uw,y_uw)
    norm_s = np.trapz([np.trapz(n, x_uw) for n in n_uw], y_uw)
    p_uw = setup_rbf_interpolator(xx0.flatten(), yy0.flatten(), (n_uw.flatten()+1e-5), logx=True, logy=True, scale=[fx,fy])

    if do_diagnostic_plot:
        import tools.plots as pl
        mcx = np.random.rand(xx.size/2)*(xx.max()-xx.min()) + xx.min()
        mcy = np.random.rand(yy.size/2)*(yy.max()-yy.min()) + yy.min()
        samples = p_tar(np.log10(mcx)/fx, np.log10(mcy)/fy)
        pl.interpolator_diagnostic(n_tar, xx.flatten(), yy.flatten(), samples, mcx, mcy)

    print "Interpolating to data."

    if not split:
        try:
            pt = p_tar(np.log10(unweighted[0])/fx, np.log10(unweighted[1])/fy)/norm
            ps = p_uw(np.log10(unweighted[0])/fx, np.log10(unweighted[1])/fy)/norm_s
            p = pt/ps *  (ps>1e-4).astype(int)
        except Exception as error:
            print "Whole array interpolation did not work. The error was", error
            split=True

    if split:
        stalle=False
        npieces = 4
        ngal=unweighted[0].size
        p = np.zeros(ngal)

        while not stalle:
            lpiece = ngal/npieces
            print "Trying %d pieces of length %3.3f M ..."%(npieces, lpiece/1e6),
            for n in xrange(npieces):
                pt = p_tar(np.log10(unweighted[0][n*lpiece : (n+1)*lpiece])/fx, np.log10(unweighted[1][n*lpiece : (n+1)*lpiece])/fy)
                ps = p_uw(np.log10(unweighted[0][n*lpiece : (n+1)*lpiece])/fx, np.log10(unweighted[1][n*lpiece : (n+1)*lpiece])/fy)
                p[n*lpiece : (n+1)*lpiece] = pt/ps * ((ps>1e-3) & (pt>1e-3)).astype(int)
            print "done."

            stalle = True
#            except:
#                print "no - will try again."
#                npieces *= 2
#                continue

    p[p<0] = 0

    import pdb ; pdb.set_trace()


    return p

def get_alpha(xdata, catalogue, nbins=5, external_calibration_col=None, apply_calibration=False, ellipticity_name="e", xdata_name="mean_hsm_psf_e%d_sky", use_weights=False, weights=None, xlim=(-1.,1.), names=["alpha","c"], binning="equal_number", silent=False, visual=False, return_vals=False):

    g1 = xdata[xdata_name%1]
    g2 = xdata[xdata_name%2]
    sel = (g1>xlim[0]) & (g1<xlim[1]) & (g2>xlim[0]) & (g2<xlim[1]) & np.isfinite(g1) & np.isfinite(g2)
    if apply_calibration and (external_calibration_col is None):
        mask_nans = np.isfinite(catalogue["m"]) & np.isfinite(catalogue["c1"]) & np.isfinite(catalogue["c2"]) 
        sel = sel & mask_nans
    g1 = g1[sel]
    g2 = g2[sel]

    e1 = catalogue["%s1"%ellipticity_name][sel]
    e2 = catalogue["%s2"%ellipticity_name][sel]

    if "weight" in catalogue.dtype.names and (use_weights):
        w = catalogue["weight"][sel]
    elif weights is not None:
        print "Will use weights provided"
        w = weights[sel]
    else:
        w = np.ones_like(e1)

    if apply_calibration:
        if external_calibration_col is None:
            m = catalogue["m"][sel]
            c1 = catalogue["c1"][sel]
            c2 = catalogue["c2"][sel]
            print "Applying calibration columns"
        else:
            m = external_calibration_col["m"][sel]
            c1 = external_calibration_col["c1"][sel]
            c2 = external_calibration_col["c2"][sel]
            print "Applying external calibration columns"
    else:
        m = np.zeros_like(e1)
        c1 = np.zeros_like(e1)
        c2 = np.zeros_like(e1)

    if  ("weight" not in catalogue.dtype.names) and (use_weights):
        print "Warning: you set use_weights=True, but there is no weights column."
        print "using unweighted values."

    if isinstance(binning,str):
        if binning=="equal_number":
            bins = find_bin_edges(g1, nbins)
        if binning=="uniform":
            bins = np.linspace(g1.min(),g1.max(), nbins+1)
    else:
        bins = binning

    x = (bins[1:]+bins[:-1])/2.0

    y11,y22,y12,y21=[],[],[],[]
    variance_y11,variance_y22,variance_y12,variance_y21=[],[],[],[]
    
    for i, bounds in enumerate(zip(bins[:-1],bins[1:])):
        lower = bounds[0]
        upper = bounds[1]
        sel1 = (g1>lower) & (g1<upper)
        sel2 = (g2>lower) & (g2<upper)

        y11.append(np.sum(w[sel1]*(e1[sel1]-c1[sel1])) / np.sum(w[sel1]*(1+m[sel1])))
        variance_y11.append( compute_im3shape_weight(e1[sel1]-g1[sel1], verbose=False) / (e1[sel1].size**0.5) )

        y22.append(np.sum(w[sel2]*(e2[sel2]-c2[sel2])) / np.sum(w[sel2]*(1+m[sel2])))
        if not (np.isfinite(np.array(y22)).all()):
            import pdb ; pdb.set_trace()
        variance_y22.append( compute_im3shape_weight(e2[sel2]-g2[sel2], verbose=False) / (e2[sel2].size**0.5) )


        y12.append(np.sum(w[sel2]*(e1[sel2]-c1[sel2])) / np.sum(w[sel2]*(1+m[sel2])))
        variance_y12.append( compute_im3shape_weight(e2[sel1]-g2[sel1], verbose=False) / (e2[sel1].size**0.5) )
       
        y21.append(np.sum(w[sel1]*(e2[sel1]-c2[sel1])) / np.sum(w[sel1]*(1+m[sel1])))
        variance_y21.append( compute_im3shape_weight(e1[sel2]-g1[sel2], verbose=False) / (e1[sel2].size**0.5) )

        

    d1 = np.array(y11)
    d2 = np.array(y22)

    # Four linear fits to do here

    # m11, c11
    p11, cov11 = optimise.curve_fit(fline, x, d1, sigma=np.array(variance_y11), p0=[0.15,1e-5])
    if not np.isfinite(cov11).all():
        p11,cov11 = stabilised_fit(x, d1, variance_y11)
    m11,c11 = p11

    # m22, c22
    p22, cov22 = optimise.curve_fit(fline, x, d2, sigma=np.array(variance_y22), p0=[0.15,1e-5])
    if not np.isfinite(cov22).all():
        p22,cov22 = stabilised_fit(x, d2, variance_y22)
    m22,c22 = p22

    # m12, c12
    p12, cov12 =  optimise.curve_fit(fline, x, y12, sigma=np.array(variance_y12), p0=[0.15,1e-5])
    if not np.isfinite(cov12).all():
        p12,cov12 = stabilised_fit(x, y12, variance_y12)
    m12,c12 = p12

    # m21, c21
    p21, cov21 =  optimise.curve_fit(fline, x, y21, sigma=np.array(variance_y21), p0=[0.15,1e-5])
    if not np.isfinite(cov21).all():
        p21,cov21 = stabilised_fit(x, y21, variance_y21)
    m21,c21 = p21

    if not silent:
        print 'alpha11=%f +- %f'%(m11,cov11[0,0]**0.5)
        print 'alpha22=%f +- %f'%(m22,cov22[0,0]**0.5)
        print 'alpha12=%f +- %f'%(m12,cov12[0,0]**0.5)
        print 'alpha21=%f +- %f'%(m21,cov21[0,0]**0.5)
        print 'c11=%f +- %f'%(c11,cov11[1,1]**0.5)
        print 'c22=%f +- %f'%(c22,cov22[1,1]**0.5)
        print 'c12=%f +- %f'%(c12,cov12[1,1]**0.5)
        print 'c21=%f +- %f'%(c21,cov21[1,1]**0.5)

    m = (m11+m22)/2
    c = (c11+c22)/2
    error_m = np.sqrt(cov11[0,0] + cov22[0,0])
    error_c = np.sqrt(cov11[1,1] + cov22[1,1])

    if visual:
        import pylab as plt
        fig, ax1 = plt.subplots()
        plt.subplots_adjust(wspace=0, hspace=0, top=0.85, bottom=0.06)
        ax1.errorbar(x,y11, variance_y11, fmt="o", color="purple")
        ax1.errorbar(x,y22, variance_y22, fmt="D", color="steelblue")
        ax1.plot(x,x*m11+c11, lw=2.5, color="purple", label=r"$\alpha_{11} = %1.4f +- %1.4f$"%(m11,cov11[0,0]**0.5))
        ax1.plot(x,x*m22+c22, lw=2.5, color="steelblue", label=r"$\alpha_{22} = %1.4f +- %1.4f$"%(m22,cov22[0,0]**0.5))
        ax1.set_xlabel("PSF Ellipticicty $e^{PSF}_{i}$")
        ax1.set_ylabel("Ellipticity $e_{i}$")
        ax1.set_xlim(xlim[0],xlim[1])
        ax1.set_ylim(-0.005,0.005)
        ax1.axhline(0,color="k", lw=2.)

        ax2 = ax1.twinx()
        ax2.set_xlim(xlim[0], xlim[1])
        plt.setp(ax2.get_yticklabels(), visible=False)
        ax2.hist(g1, alpha=0.1, bins=45, color="purple")
        ax2.hist(g2, alpha=0.1, bins=45, color="steelblue")
        ax1.legend(loc="upper left")

        plt.tight_layout()


    biases_dict = {"alpha":(m,error_m), "c":(c,error_c), "alpha11":(m11,cov11[0,0]**0.5), "c11":(c11,cov11[1,1]**0.5), "alpha22":(m22,cov22[0,0]**0.5), "c22":(c22,cov22[1,1]**0.5), "alpha12":(m12,cov12[0,0]**0.5), "c12":(c12,cov12[1,1]**0.5), "alpha21":(m21,cov21[0,0]**0.5), "c21":(c21,cov21[1,1]**0.5)}

    # Decide how to package the results
    # Either as a dictionary, including the fit covariances
    # Or as a single value (probably for when using bootstrap errors)
    out={}

    if not isinstance(names, str):
        for name in names:
            out[name] = biases_dict[name]
    else:
        out = biases_dict[names][0]

    if return_vals:
        out["e1"] = y11
        out["de1"] = variance_y11
        out["e2"] = y22
        out["de2"] = variance_y22
        out["bins"] = x

    return out

def fline(x, m, c):
    return m*x + c

def find_distance_self(array, name="cosmos_ident", check_array=[], verbose=True):
    import scipy.spatial as sps
    R = [] #np.zeros_like(array[name]).astype(float)
    ids = np.unique(array[name])

    ind=[]


    print "Found %d unique ids."%ids.size
    for n, i in enumerate(ids):
        if len(check_array)!=0:
            if array["DES_id"][array[name]==i][0] in check_array["coadd_objects_id"]:
                continue
        #select = np.argwhere(array[name]==i).T[0]
        xy = np.vstack((array["ra"][array[name]==i], array["dec"][array[name]==i]))

        if xy.shape[1]<2:
            result = [9000.]
        else:
            tree = sps.KDTree(xy.T)
            result = tree.query(xy.T, k=2)[0].T[1]

        with open("self_distance.txt", "a") as f:
            out = np.vstack((array["DES_id"][array[name]==i], result))
            np.savetxt(f, out.T)


        R+=list(result)
        ind+=list(array[array[name]==i])

       # R[select] = result

        if verbose:
            print n, i

    R = np.array(R)*60*60/0.27
    ind = np.array(ind).astype(int)
#
    dt = [(name, int),("self_distance", float)]
    out = np.zeros(R.size, dtype=dt)
    out[name] = array[name]
    out["self_distance"] = R
#
    file = fi.FITS("cosmos_self_distance.fits", "rw")
    file.write(out, clobber=True)
    file.close()

    return R

def extract_cosmos_column(cosmos, target, column_name, outfile=None, start_point=0, return_vals=True):
    ids = np.unique(target["cosmos_ident"])
    column = np.zeros_like(target["cosmos_ident"]).astype(float)

    if outfile is None:
        outfile = "%s.txt"%column_name

    if os.path.exists(outfile):
        reference = open(outfile).read()
    else:
        reference=[]

    print "Writing to disc %s"%outfile
    f = open(outfile, "wa")

    print "Found %d unique ids."%ids.size

    print "Will start from element %d"%start_point

    for n, i in enumerate(target["cosmos_ident"][start_point:]):
        if "%d"%target["DES_id"][n] in reference:
            print "%d %d -- already done "%(n, i)
            continue 
        val = cosmos[cosmos["ident"]==i][0][column_name]
        if return_vals: 
            column[(target["cosmos_ident"]==i)] = val

        f.write("%d %f \n"%(target["DES_id"][n],val))

        print n, i
    f.close()

    if return_vals:
        print "packaging array"

        out = np.zeros(len(column), dtype=(column[0].dtype, "%s"%column_name))
        out[column_name] = column

        return out
        




def get_correct_boxsize(cosmos, results, truth):
    f = np.sqrt(-2*np.log(0.5))
    for i, row in enumerate(results):
        print i
        c = cosmos[cosmos["ident"]==data["cosmos_ident"]]
        rp = f * results["mean_hsm_psf_sigma"]/3

def identity(inp):
    return inp

def setup_rbf_interpolator(x, y, z, logx=True, logy=True, scale=[]):

    if logx:
        xfunction = np.log10
    else:
        xfunction = di.identity
    if logy:
        yfunction = np.log10
    else:
        yfunction = di.identity
        
    # Set up the RBF interpolation
    if len(scale)>0:
        fx = xfunction(scale[0])
        fy = yfunction(scale[1])
    else:
        fx = fy = 1

    return Rbf(xfunction(x)/fx, yfunction(y)/fy, z, smooth=0, function="multiquadric")
        


#
#f=np.sqrt(-2*np.log(0.5))                    # HLR - sigma conversion
#fac=2.35                                     # sigma - FWHM conversion
#R=hoopoe.truth["hlr"]/f                      #Gaussian sigma galaxy, arcsec
#rp=hoopoe.res["mean_hsm_psf_sigma"]*0.27/3   # sigma PSF, factor of 3 for the upsampling, arcsec
#
#sig=np.sqrt(R*R + rp*rp)/0.27                # sigma convolved, pixels                                
#eps = (1-q0)
#b = 2*5* sig * (1+eps)
#
#boxsize=np.zeros_like(q0)
#for b0,b1 in zip(allowed_boxsizes[:-1], allowed_boxsizes[1:]):
#    select = (b>b0) & (b<b1)
#    boxsize[select] = b1
#    print b1


relevant_parameters = ['ra_as', 'dec_as', 'e1', 'e2', 'radius', 'radius_ratio', 'bulge_a', 'disc_a', 'coadd_objects_id', 'time', 'bulge_flux', 'disc_flux', 'flux_ratio', 'snr', 'old_snr', 'min_residuals', 'max_residuals', 'model_min', 'model_max', 'likelihood', 'levmar_start_error', 'levmar_end_error', 'levmar_resid_grad', 'levmar_vector_diff', 'levmar_error_diff', 'levmar_comp_grad', 'levmar_iterations', 'levmar_reason', 'levmar_like_evals', 'levmar_grad_evals', 'levmar_sys_evals', 'mean_flux', 'n_exposure', 'stamp_size', 'mean_rgpp_rp', 'mean_psf_e1_sky', 'mean_psf_e2_sky', 'fails_psf_e2_sky', 'mean_psf_fwhm', 'mean_unmasked_flux_frac', 'fails_unmasked_flux_frac', 'mean_model_edge_mu', 'mean_model_edge_sigma', 'mean_edge_mu', 'fails_edge_mu', 'mean_edge_sigma', 'mean_hsm_psf_e1_sky', 'mean_hsm_psf_e2_sky', 'mean_hsm_psf_sigma', 'mean_hsm_psf_rho4', 'mean_mask_fraction',  'round_snr',  'round_snr_mw', 'ra', 'dec', 'chi2_pixel', 'mag_auto_g', 'mag_auto_r', 'mag_auto_i', 'mag_auto_z', 'desdm_zp']
    
i3sdt=np.dtype([('ra_as', '>f8'), ('dec_as', '>f8'), ('e1', '>f8'), ('e2', '>f8'), ('radius', '>f8'), ('radius_ratio', '>f8'), ('bulge_A', '>f8'), ('disc_A', '>f8'), ('bulge_index', '>f8'), ('disc_index', '>f8'), ('identifier', 'i8'), ('time', '>f8'), ('bulge_flux', '>f8'), ('disc_flux', '>f8'), ('flux_ratio', '>f8'), ('snr', '>f8'), ('old_snr', '>f8'), ('min_residuals', '>f8'), ('max_residuals', '>f8'), ('model_min', '>f8'), ('model_max', '>f8'), ('likelihood', '>f8'), ('levmar_start_error', '>f8'), ('levmar_end_error', '>f8'), ('levmar_resid_grad', '>f8'), ('levmar_vector_diff', '>f8'), ('levmar_error_diff', '>f8'), ('levmar_comp_grad', '>f8'), ('levmar_iterations', '>f8'), ('levmar_reason', '>f8'), ('levmar_like_evals', '>f8'), ('levmar_grad_evals', '>f8'), ('levmar_sys_evals', '>f8'), ('mean_flux', '>f8'), ('number_varied_params', '>f8'), ('covmat_0_0', '>f8'), ('covmat_0_1', '>f8'), ('covmat_0_2', '>f8'), ('covmat_0_3', '>f8'), ('covmat_0_4', '>f8'), ('covmat_0_5', '>f8'), ('covmat_1_0', '>f8'), ('covmat_1_1', '>f8'), ('covmat_1_2', '>f8'), ('covmat_1_3', '>f8'), ('covmat_1_4', '>f8'), ('covmat_1_5', '>f8'), ('covmat_2_0', '>f8'), ('covmat_2_1', '>f8'), ('covmat_2_2', '>f8'), ('covmat_2_3', '>f8'), ('covmat_2_4', '>f8'), ('covmat_2_5', '>f8'), ('covmat_3_0', '>f8'), ('covmat_3_1', '>f8'), ('covmat_3_2', '>f8'), ('covmat_3_3', '>f8'), ('covmat_3_4', '>f8'), ('covmat_3_5', '>f8'), ('covmat_4_0', '>f8'), ('covmat_4_1', '>f8'), ('covmat_4_2', '>f8'), ('covmat_4_3', '>f8'), ('covmat_4_4', '>f8'), ('covmat_4_5', '>f8'), ('covmat_5_0', '>f8'), ('covmat_5_1', '>f8'), ('covmat_5_2', '>f8'), ('covmat_5_3', '>f8'), ('covmat_5_4', '>f8'), ('covmat_5_5', '>f8'), ('n_exposure', '>f8'), ('nparam_varied', '>f8'), ('stamp_size', '>f8'), ('mean_rgpp_rp', '>f8'), ('fails_rgpp_rp', 'i8'), ('mean_psf_e1_sky', '>f8'), ('fails_psf_e1_sky', 'i8'), ('mean_psf_e2_sky', '>f8'), ('fails_psf_e2_sky', 'i8'), ('mean_psf_fwhm', '>f8'), ('fails_psf_fwhm', '>f8'), ('mean_unmasked_flux_frac', '>f8'), ('fails_unmasked_flux_frac', '>f8'), ('mean_model_edge_mu', '>f8'), ('fails_model_edge_mu', 'i8'), ('mean_model_edge_sigma', '>f8'), ('fails_model_edge_sigma', 'i8'), ('mean_edge_mu', '>f8'), ('fails_edge_mu', 'i8'), ('mean_edge_sigma', '>f8'), ('fails_edge_sigma', '>f8'), ('mean_hsm_psf_e1_sky', '>f8'), ('fails_hsm_psf_e1_sky', 'i8'), ('mean_hsm_psf_e2_sky', '>f8'), ('fails_hsm_psf_e2_sky', 'i8'), ('mean_hsm_psf_sigma', '>f8'), ('fails_hsm_psf_sigma', '>f8'), ('mean_hsm_psf_rho4', '>f8'), ('fails_hsm_psf_rho4', 'i8'), ('mean_mask_fraction', '>f8'), ('fails_mask_fraction', 'i8'), ('round_snr', '>f8'), ('fails_round_snr', 'i8'), ('round_snr_mw', '>f8'), ('fails_round_snr_mw', 'i8'), ('bands', 'S5'), ('tilename', 'S20'), ('ra', '>f8'), ('row_id', '>f8'), ('dec', '>f8'), ('id', 'i8')])

ppi3sdt=np.dtype([('ra_as', '>f8'), ('dec_as', '>f8'), ('e1', '>f8'), ('e2', '>f8'), ('radius', '>f8'), ('radius_ratio', '>f8'), ('bulge_a', '>f8'), ('disc_a', '>f8'), ('bulge_index', '>f8'), ('disc_index', '>f8'), ('coadd_objects_id', '>i8'), ('time', '>f8'), ('bulge_flux', '>f8'), ('disc_flux', '>f8'), ('flux_ratio', '>f8'), ('snr', '>f8'), ('old_snr', '>f8'), ('min_residuals', '>f8'), ('max_residuals', '>f8'), ('model_min', '>f8'), ('model_max', '>f8'), ('likelihood', '>f8'), ('levmar_start_error', '>f8'), ('levmar_end_error', '>f8'), ('levmar_resid_grad', '>f8'), ('levmar_vector_diff', '>f8'), ('levmar_error_diff', '>f8'), ('levmar_comp_grad', '>f8'), ('levmar_iterations', '>i8'), ('levmar_reason', '>i8'), ('levmar_like_evals', '>i8'), ('levmar_grad_evals', '>i8'), ('levmar_sys_evals', '>i8'), ('mean_flux', '>f8'), ('number_varied_params', '>i8'), ('covmat_0_0', '>f8'), ('covmat_0_1', '>f8'), ('covmat_0_2', '>f8'), ('covmat_0_3', '>f8'), ('covmat_0_4', '>f8'), ('covmat_0_5', '>f8'), ('covmat_1_0', '>f8'), ('covmat_1_1', '>f8'), ('covmat_1_2', '>f8'), ('covmat_1_3', '>f8'), ('covmat_1_4', '>f8'), ('covmat_1_5', '>f8'), ('covmat_2_0', '>f8'), ('covmat_2_1', '>f8'), ('covmat_2_2', '>f8'), ('covmat_2_3', '>f8'), ('covmat_2_4', '>f8'), ('covmat_2_5', '>f8'), ('covmat_3_0', '>f8'), ('covmat_3_1', '>f8'), ('covmat_3_2', '>f8'), ('covmat_3_3', '>f8'), ('covmat_3_4', '>f8'), ('covmat_3_5', '>f8'), ('covmat_4_0', '>f8'), ('covmat_4_1', '>f8'), ('covmat_4_2', '>f8'), ('covmat_4_3', '>f8'), ('covmat_4_4', '>f8'), ('covmat_4_5', '>f8'), ('covmat_5_0', '>f8'), ('covmat_5_1', '>f8'), ('covmat_5_2', '>f8'), ('covmat_5_3', '>f8'), ('covmat_5_4', '>f8'), ('covmat_5_5', '>f8'), ('n_exposure', '>i8'), ('nparam_varied', '>i8'), ('stamp_size', '>i8'), ('mean_rgpp_rp', '>f8'), ('fails_rgpp_rp', '>i8'), ('mean_psf_e1_sky', '>f8'), ('fails_psf_e1_sky', '>i8'), ('mean_psf_e2_sky', '>f8'), ('fails_psf_e2_sky', '>i8'), ('mean_psf_fwhm', '>f8'), ('fails_psf_fwhm', '>i8'), ('mean_unmasked_flux_frac', '>f8'), ('fails_unmasked_flux_frac', '>i8'), ('mean_model_edge_mu', '>f8'), ('fails_model_edge_mu', '>i8'), ('mean_model_edge_sigma', '>f8'), ('fails_model_edge_sigma', '>i8'), ('mean_edge_mu', '>f8'), ('fails_edge_mu', '>i8'), ('mean_edge_sigma', '>f8'), ('fails_edge_sigma', '>i8'), ('mean_hsm_psf_e1_sky', '>f8'), ('fails_hsm_psf_e1_sky', '>i8'), ('mean_hsm_psf_e2_sky', '>f8'), ('fails_hsm_psf_e2_sky', '>i8'), ('mean_hsm_psf_sigma', '>f8'), ('fails_hsm_psf_sigma', '>i8'), ('mean_hsm_psf_rho4', '>f8'), ('fails_hsm_psf_rho4', '>i8'), ('mean_mask_fraction', '>f8'), ('fails_mask_fraction', '>i8'), ('round_snr', '>f8'), ('fails_round_snr', '>i8'), ('round_snr_mw', '>f8'), ('fails_round_snr_mw', '>i8'), ('bands', 'S5'), ('tilename', 'S12'), ('ra', '>f8'), ('row_id', '>i8'), ('dec', '>f8'), ('id', '>i8'), ('chi2_pixel', '>f8'), ('error_flag', '>i8'), ('info_flag', '>i8'), ('is_bulge', '>i8')])
