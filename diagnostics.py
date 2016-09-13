import numpy as np
import astropy.table
import matplotlib.pyplot as plt
import glob, os, pyfits, pdb 
import fitsio as fio

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





def load_results(res_path='None', keyword='fornax', format='txt', postprocessed=True, return_filelist=False, cols=None, ntot=-1, apply_infocuts=False):
    files=[]
    ind=[]
    Nf=0
    if apply_infocuts:
        print "Applying info flag cuts on loading."
    if format=='txt':
        extra='.main'
    else:
        extra=''

    if postprocessed:
        dt = ppi3sdt
    else:
        dt = i3sdt
    
    if cols:
        newdt = [(n, str(dt[i])) for i, n in enumerate(dt.names) if n in cols ]
        dt = np.dtype(newdt)

  # Generate empty array to fill
    if apply_infocuts:
        buff=20000
    else:
        buff=4100
    if ntot!=-1:
        res= np.empty(ntot*buff, dtype=dt)
    else:
        res= np.empty(3000*buff, dtype=dt)

    files = glob.glob('*%s*%s.%s'%(keyword,extra,format))
    if res_path!='None':
        files = glob.glob(os.path.join(res_path,'*%s*%s.%s'%(keyword,extra,format)))
    if len(files)==0:
        files = glob.glob('meds/*%s*%s.%s'%(keyword,extra,format))
    if len(files)==0:
        files = glob.glob('*%s*%s.%s'%(keyword,extra,format))

    if ntot!=-1:
        files=files[:ntot]
            
    for f in files:
        if format.lower()=="fits":
            fits = fio.FITS(f)
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

        i0 = np.argwhere(res["id"]==0)[0,0]
        i1 = i0 + len(dat)
        res[i0:i1] = dat
        ind+=[(i0,i1)]
        if res[i1-1]!=dat[-1]:
            print "ERROR: Data overflow. %d %d"%(i1, len(res))
        print Nf, f
        Nf+=1

    res = res[res["id"]!=0]

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

def load_truth(truth_path=None, keyword='DES', match=None, cols=None, ind=None, res=None):
    files=[]
    if truth_path!=None:
        files = glob.glob(os.path.join(truth_path,'*%s*-truth*.fits*'%keyword))
        if len(files)==0:
            files = glob.glob(truth_path)
    if len(files)==0:
        files = glob.glob('../truth/*%s*-truth*.fits*'%keyword)

    if match:
        list_res=[os.path.basename(m)[:12] for m in match]
        print "matching to tilelist ", list_res
            
        filelist=np.empty(len(list_res),dtype="S100")
        for f in files:
            if os.path.basename(f)[:12] in list_res:
                i=np.argwhere(np.array(list_res)==os.path.basename(f)[:12])[0,0]
                filelist[i]=f


        #filelist = np.array([np.array([ f, np.argwhere(np.array(list_res)==os.path.basename(f)[:12])[0,0] ]) for f in files if os.path.basename(f)[:12] in list_res])
    else:
        filelist=files
 
    dt = fio.read(filelist[0]).dtype
    truth = np.empty(len(filelist)*10000, dtype=dt)
    for i, f in enumerate(filelist):
        fits = fio.FITS(f) 
        if cols:
            dat = fits[1].read(columns=cols)
        else:
            dat = fits[1].read()
        fits.close()
            
        if ind!=None and res!=None:
            selres=res[ind[i][0]:ind[i][1]]
            r,dat=match_results(selres,dat)
        i0 = np.argwhere(truth['DES_id']==0)[0,0]
        i1 = i0 + len(dat)
        truth[i0:i1] = dat
    #truth = np.concatenate((np.array(truth), np.array(astropy.table.Table.read(f, format="fits"))))
    print i+1, f, "%d/%d"%(i1, len(filelist)*10000)
                
    truth = truth[truth["DES_id"]!=0]

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

def match_results(res,tr, name1="DES_id", name2="id", unique=True):

    if unique:
        un,ind=np.unique(res[name2], return_index=True)
        res=res[ind]
        un,ind=np.unique(tr[name1], return_index=True)
        tr=tr[ind]

    tr_id = tr[name1]
    res_id = res[name2]

    tr_order = np.argsort(tr_id)
    res_order = np.argsort(res_id)

    tr = tr[tr_order]
    res = res[res_order]

    if (res.shape==tr.shape) and (np.all(res[name2]==tr[name1])):
        print 'Matched %d objects'%len(tr[name1])
    else:
        sel1 = np.in1d(tr[name1], res[name2])

        tr = tr[sel1]

        sel2 = np.in1d(res[name2], tr[name1])

        res = res[sel2]

    return res,tr 

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

def bootstrap_error(nsubsamples, full_cat, operation, additional_args=None, additional_argvals=None):
    """Generic function which takes an array and splits into nsubsamples parts. Then apply an operation
       to each, and find the rms deviation about the mean in the resulting quantity."""
    resampled = []
    resample_length=len(full_cat)/nsubsamples
    bootstrap_edges=(np.arange(nsubsamples+1)*resample_length).astype(int)
    for i in xrange(nsubsamples-1):
        b_low = bootstrap_edges[i]
        b_high = bootstrap_edges[i+1]
        print "bootstrap subsample %d (%d-%d)"%(i+1, b_high, b_low)
        if additional_args is None:
            derived_quantity = operation( full_cat[b_low:b_high])
        else:
            comm = "derived_quantity = operation(full_cat[b_low:b_high]"
            for name,val in zip(additional_args, additional_argvals):
                if isinstance(val,str):
                    comm += ", %s='%s'"%(name,val)
                else:
                    comm += ", %s=%s"%(name,val)
            comm+=")"

            exec comm
            print derived_quantity

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


relevant_parameters = ['ra_as', 'dec_as', 'e1', 'e2', 'radius', 'radius_ratio', 'bulge_a', 'disc_a', 'coadd_objects_id', 'time', 'bulge_flux', 'disc_flux', 'flux_ratio', 'snr', 'old_snr', 'min_residuals', 'max_residuals', 'model_min', 'model_max', 'likelihood', 'levmar_start_error', 'levmar_end_error', 'levmar_resid_grad', 'levmar_vector_diff', 'levmar_error_diff', 'levmar_comp_grad', 'levmar_iterations', 'levmar_reason', 'levmar_like_evals', 'levmar_grad_evals', 'levmar_sys_evals', 'mean_flux', 'n_exposure', 'stamp_size', 'mean_rgpp_rp', 'mean_psf_e1_sky', 'mean_psf_e2_sky', 'fails_psf_e2_sky', 'mean_psf_fwhm', 'mean_unmasked_flux_frac', 'fails_unmasked_flux_frac', 'mean_model_edge_mu', 'mean_model_edge_sigma', 'mean_edge_mu', 'fails_edge_mu', 'mean_edge_sigma', 'mean_hsm_psf_e1_sky', 'mean_hsm_psf_e2_sky', 'mean_hsm_psf_sigma', 'mean_hsm_psf_rho4', 'mean_mask_fraction',  'round_snr',  'round_snr_mw', 'ra', 'dec', 'chi2_pixel', 'mag_auto_g', 'mag_auto_r', 'mag_auto_i', 'mag_auto_z', 'desdm_zp']
    
i3sdt=np.dtype([('ra_as', '>f8'), ('dec_as', '>f8'), ('e1', '>f8'), ('e2', '>f8'), ('radius', '>f8'), ('radius_ratio', '>f8'), ('bulge_A', '>f8'), ('disc_A', '>f8'), ('bulge_index', '>f8'), ('disc_index', '>f8'), ('identifier', 'i8'), ('time', '>f8'), ('bulge_flux', '>f8'), ('disc_flux', '>f8'), ('flux_ratio', '>f8'), ('snr', '>f8'), ('old_snr', '>f8'), ('min_residuals', '>f8'), ('max_residuals', '>f8'), ('model_min', '>f8'), ('model_max', '>f8'), ('likelihood', '>f8'), ('levmar_start_error', '>f8'), ('levmar_end_error', '>f8'), ('levmar_resid_grad', '>f8'), ('levmar_vector_diff', '>f8'), ('levmar_error_diff', '>f8'), ('levmar_comp_grad', '>f8'), ('levmar_iterations', '>f8'), ('levmar_reason', '>f8'), ('levmar_like_evals', '>f8'), ('levmar_grad_evals', '>f8'), ('levmar_sys_evals', '>f8'), ('mean_flux', '>f8'), ('number_varied_params', '>f8'), ('covmat_0_0', '>f8'), ('covmat_0_1', '>f8'), ('covmat_0_2', '>f8'), ('covmat_0_3', '>f8'), ('covmat_0_4', '>f8'), ('covmat_0_5', '>f8'), ('covmat_1_0', '>f8'), ('covmat_1_1', '>f8'), ('covmat_1_2', '>f8'), ('covmat_1_3', '>f8'), ('covmat_1_4', '>f8'), ('covmat_1_5', '>f8'), ('covmat_2_0', '>f8'), ('covmat_2_1', '>f8'), ('covmat_2_2', '>f8'), ('covmat_2_3', '>f8'), ('covmat_2_4', '>f8'), ('covmat_2_5', '>f8'), ('covmat_3_0', '>f8'), ('covmat_3_1', '>f8'), ('covmat_3_2', '>f8'), ('covmat_3_3', '>f8'), ('covmat_3_4', '>f8'), ('covmat_3_5', '>f8'), ('covmat_4_0', '>f8'), ('covmat_4_1', '>f8'), ('covmat_4_2', '>f8'), ('covmat_4_3', '>f8'), ('covmat_4_4', '>f8'), ('covmat_4_5', '>f8'), ('covmat_5_0', '>f8'), ('covmat_5_1', '>f8'), ('covmat_5_2', '>f8'), ('covmat_5_3', '>f8'), ('covmat_5_4', '>f8'), ('covmat_5_5', '>f8'), ('n_exposure', '>f8'), ('nparam_varied', '>f8'), ('stamp_size', '>f8'), ('mean_rgpp_rp', '>f8'), ('fails_rgpp_rp', 'i8'), ('mean_psf_e1_sky', '>f8'), ('fails_psf_e1_sky', 'i8'), ('mean_psf_e2_sky', '>f8'), ('fails_psf_e2_sky', 'i8'), ('mean_psf_fwhm', '>f8'), ('fails_psf_fwhm', '>f8'), ('mean_unmasked_flux_frac', '>f8'), ('fails_unmasked_flux_frac', '>f8'), ('mean_model_edge_mu', '>f8'), ('fails_model_edge_mu', 'i8'), ('mean_model_edge_sigma', '>f8'), ('fails_model_edge_sigma', 'i8'), ('mean_edge_mu', '>f8'), ('fails_edge_mu', 'i8'), ('mean_edge_sigma', '>f8'), ('fails_edge_sigma', '>f8'), ('mean_hsm_psf_e1_sky', '>f8'), ('fails_hsm_psf_e1_sky', 'i8'), ('mean_hsm_psf_e2_sky', '>f8'), ('fails_hsm_psf_e2_sky', 'i8'), ('mean_hsm_psf_sigma', '>f8'), ('fails_hsm_psf_sigma', '>f8'), ('mean_hsm_psf_rho4', '>f8'), ('fails_hsm_psf_rho4', 'i8'), ('mean_mask_fraction', '>f8'), ('fails_mask_fraction', 'i8'), ('round_snr', '>f8'), ('fails_round_snr', 'i8'), ('round_snr_mw', '>f8'), ('fails_round_snr_mw', 'i8'), ('bands', 'S5'), ('tilename', 'S20'), ('ra', '>f8'), ('row_id', '>f8'), ('dec', '>f8'), ('id', 'i8')])

ppi3sdt=np.dtype([('ra_as', '>f8'), ('dec_as', '>f8'), ('e1', '>f8'), ('e2', '>f8'), ('radius', '>f8'), ('radius_ratio', '>f8'), ('bulge_a', '>f8'), ('disc_a', '>f8'), ('bulge_index', '>f8'), ('disc_index', '>f8'), ('coadd_objects_id', '>i8'), ('time', '>f8'), ('bulge_flux', '>f8'), ('disc_flux', '>f8'), ('flux_ratio', '>f8'), ('snr', '>f8'), ('old_snr', '>f8'), ('min_residuals', '>f8'), ('max_residuals', '>f8'), ('model_min', '>f8'), ('model_max', '>f8'), ('likelihood', '>f8'), ('levmar_start_error', '>f8'), ('levmar_end_error', '>f8'), ('levmar_resid_grad', '>f8'), ('levmar_vector_diff', '>f8'), ('levmar_error_diff', '>f8'), ('levmar_comp_grad', '>f8'), ('levmar_iterations', '>i8'), ('levmar_reason', '>i8'), ('levmar_like_evals', '>i8'), ('levmar_grad_evals', '>i8'), ('levmar_sys_evals', '>i8'), ('mean_flux', '>f8'), ('number_varied_params', '>i8'), ('covmat_0_0', '>f8'), ('covmat_0_1', '>f8'), ('covmat_0_2', '>f8'), ('covmat_0_3', '>f8'), ('covmat_0_4', '>f8'), ('covmat_0_5', '>f8'), ('covmat_1_0', '>f8'), ('covmat_1_1', '>f8'), ('covmat_1_2', '>f8'), ('covmat_1_3', '>f8'), ('covmat_1_4', '>f8'), ('covmat_1_5', '>f8'), ('covmat_2_0', '>f8'), ('covmat_2_1', '>f8'), ('covmat_2_2', '>f8'), ('covmat_2_3', '>f8'), ('covmat_2_4', '>f8'), ('covmat_2_5', '>f8'), ('covmat_3_0', '>f8'), ('covmat_3_1', '>f8'), ('covmat_3_2', '>f8'), ('covmat_3_3', '>f8'), ('covmat_3_4', '>f8'), ('covmat_3_5', '>f8'), ('covmat_4_0', '>f8'), ('covmat_4_1', '>f8'), ('covmat_4_2', '>f8'), ('covmat_4_3', '>f8'), ('covmat_4_4', '>f8'), ('covmat_4_5', '>f8'), ('covmat_5_0', '>f8'), ('covmat_5_1', '>f8'), ('covmat_5_2', '>f8'), ('covmat_5_3', '>f8'), ('covmat_5_4', '>f8'), ('covmat_5_5', '>f8'), ('n_exposure', '>i8'), ('nparam_varied', '>i8'), ('stamp_size', '>i8'), ('mean_rgpp_rp', '>f8'), ('fails_rgpp_rp', '>i8'), ('mean_psf_e1_sky', '>f8'), ('fails_psf_e1_sky', '>i8'), ('mean_psf_e2_sky', '>f8'), ('fails_psf_e2_sky', '>i8'), ('mean_psf_fwhm', '>f8'), ('fails_psf_fwhm', '>i8'), ('mean_unmasked_flux_frac', '>f8'), ('fails_unmasked_flux_frac', '>i8'), ('mean_model_edge_mu', '>f8'), ('fails_model_edge_mu', '>i8'), ('mean_model_edge_sigma', '>f8'), ('fails_model_edge_sigma', '>i8'), ('mean_edge_mu', '>f8'), ('fails_edge_mu', '>i8'), ('mean_edge_sigma', '>f8'), ('fails_edge_sigma', '>i8'), ('mean_hsm_psf_e1_sky', '>f8'), ('fails_hsm_psf_e1_sky', '>i8'), ('mean_hsm_psf_e2_sky', '>f8'), ('fails_hsm_psf_e2_sky', '>i8'), ('mean_hsm_psf_sigma', '>f8'), ('fails_hsm_psf_sigma', '>i8'), ('mean_hsm_psf_rho4', '>f8'), ('fails_hsm_psf_rho4', '>i8'), ('mean_mask_fraction', '>f8'), ('fails_mask_fraction', '>i8'), ('round_snr', '>f8'), ('fails_round_snr', '>i8'), ('round_snr_mw', '>f8'), ('fails_round_snr_mw', '>i8'), ('bands', 'S5'), ('tilename', 'S12'), ('ra', '>f8'), ('row_id', '>i8'), ('dec', '>f8'), ('id', '>i8'), ('chi2_pixel', '>f8'), ('error_flag', '>i8'), ('info_flag', '>i8'), ('is_bulge', '>i8')])
