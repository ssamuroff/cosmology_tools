import numpy as np
import fitsio as fi
import os,yaml,glob
import argparse
from string import Template as Tmp


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--nofz', type=str, action='store')
parser.add_argument('--cuts', type=str, default='none', action='store')
parser.add_argument('--outdir', type=str, action='store')

args = parser.parse_args()

comm = "python -m tools.des.pz.cosmos_resampled --cuts %s --output %s"%(args.cuts,args.outdir)
print comm
os.system(comm)


#dat = fi.FITS(args.nofz)
#nz = [dat["nz_source"]["bin%d"%i][:] for i in [1,2,3,4]]
#z = dat["nz_source"]["z_mid"][:]

#get_data_means(args)

#data_zmeans = [np.trapz(n*z,z)/np.trapz(n,z) for n in nz ]
#print "Measured BPZ bin means from data:"
#print data_zmeans

cosmos_zmeans = np.loadtxt("%s/cosmos-zmean-bins-mof_iprior.txt"%args.outdir)[1:]
print "Measured BPZ bin means from resampled COSMOS:"

def get_data_means(args):
    f = "/share/des/disc7/samuroff/des/mcal-y1a1-combined-riz-blind-v4-matched.fits"
    print "loading shapes file ", f
    rphot = fi.read(f) 
    print "Loading metacal photo-z file"
    rmcalbpz=fi.read("/share/des/disc4/samuroff/ia/pz_calibration/cats/mcal-y1a1-combined-griz-blind-v3-matched_BPZ.fits") # BPZ run on MCAL photometry
    print "Loading MOF photo-z file"
    rmofbpz=fi.read("/share/des/disc7/samuroff/des/y1a1-gold-mof-badregion_BPZ.fits")

    flags = (rphot["flags_select"]==0) # & (rphot["flags_select_1p"]==0) & (rphot["flags_select_2p"]==0) & (rphot["flags_select_1m"]==0) & (rphot["flags_select_2m"]==0)
    rmcalbpz = rmcalbpz[flags]
    rmofbpz = rmofbpz[flags]
    rphot = rphot[flags]

    print(np.max(rphot['coadd_objects_id']-rmcalbpz['coadd_objects_id'])) # check IDs are the same

    if args.cuts!='none':
        cuts = yaml.load(open(args.cuts))
    else:
       cuts = {}
    mcal_mask = np.ones(rmcalbpz.size).astype(bool)
    mof_mask = np.ones(rmofbpz.size).astype(bool)
    for col in cuts.keys():
        print "Cutting %s"%col
        lower,upper,dt = cuts[col].split()
        lower=float(lower) ; upper=float(upper)
        mcal_mask = mcal_mask & (rmcalbpz[col.lower()]>lower) & (rmcalbpz[col.lower()]<upper)
        mof_mask = mof_mask & (rmofbpz[col.lower()]>lower) & (rmofbpz[col.lower()]<upper)

        # assign to bins in MCAL MEAN_Z
        # 0=all: 0.2<mean_z<1.3
        # 1=0.2..0.43
        # 2=0.43..0.63
        # 3=0.63..0.90
        # 4=0.90..1.30
    lmask0,=np.where((rmcalbpz['mean_z']>0.20)&(rmcalbpz['mean_z']<1.30) & mcal_mask)
    lmask1,=np.where((rmcalbpz['mean_z']>0.20)&(rmcalbpz['mean_z']<0.43) & mcal_mask)
    lmask2,=np.where((rmcalbpz['mean_z']>0.43)&(rmcalbpz['mean_z']<0.63) & mcal_mask)
    lmask3,=np.where((rmcalbpz['mean_z']>0.63)&(rmcalbpz['mean_z']<0.90) & mcal_mask)
    lmask4,=np.where((rmcalbpz['mean_z']>0.90)&(rmcalbpz['mean_z']<1.30) & mcal_mask)

    mask0,=np.where((rmcalbpz['mean_z']>0.20)&(rmcalbpz['mean_z']<1.30))
    mask1,=np.where((rmcalbpz['mean_z']>0.20)&(rmcalbpz['mean_z']<0.43))
    mask2,=np.where((rmcalbpz['mean_z']>0.43)&(rmcalbpz['mean_z']<0.63))
    mask3,=np.where((rmcalbpz['mean_z']>0.63)&(rmcalbpz['mean_z']<0.90))
    mask4,=np.where((rmcalbpz['mean_z']>0.90)&(rmcalbpz['mean_z']<1.30))

    def wavg(x,wt):
        return sum(x*wt)/sum(wt)

    print("all <z>=",  wavg(rmofbpz['z_mc'][lmask0],rphot['R11'][lmask0]+rphot['R22'][lmask0]))
    print("bin 1 <z>=",wavg(rmofbpz['z_mc'][lmask1],rphot['R11'][lmask1]+rphot['R22'][lmask1]))
    print("bin 2 <z>=",wavg(rmofbpz['z_mc'][lmask2],rphot['R11'][lmask2]+rphot['R22'][lmask2]))
    print("bin 3 <z>=",wavg(rmofbpz['z_mc'][lmask3],rphot['R11'][lmask3]+rphot['R22'][lmask3]))
    print("bin 4 <z>=",wavg(rmofbpz['z_mc'][lmask4],rphot['R11'][lmask4]+rphot['R22'][lmask4]))
    print("these numbers go to https://docs.google.com/document/d/1Bo_zMI1S2F-Han7KkxAS-tDHErKJq0nt9YPu9Vl6EKc")
    print("\"COSMOS <z>, metacal\"")

    zavg = np.array([wavg(rmofbpz['z_mc'][lm],rphot['R11'][lm]+rphot['R22'][lm]) for lm in [lmask0,lmask1,lmask2,lmask3,lmask4]])
    frac = np.array([rmofbpz['z_mc'][lm[0]].size * 1.0 / rmofbpz['z_mc'][lm[1]].size for lm in [(lmask0,mask0),(lmask1,mask1),(lmask2,mask2),(lmask3,mask3),(lmask4,mask4)]])
    return zavg, frac

def export_priors(centres, widths, filename):
    template = Tmp(open("/home/samuroff/local/python/lib/python2.7/site-packages/tools/des/pz/mcal_priors_template.ini").read())
    priors_info = template.substitute(Z1=centres[0],Z2=centres[1],Z3=centres[2],Z4=centres[3],DZ1=widths[0],DZ2=widths[1],DZ3=widths[2],DZ4=widths[3])

    outfile = open(filename, "wa")
    outfile.write(priors_info)
    outfile.close()
    print "Wrote new priors to %s"%filename

data_zmeans, frac = get_data_means(args)
dz = cosmos_zmeans - data_zmeans[1:]
frac = frac[1:] 

comm = "python -m tools.des.pz.get_prior_widths --cuts %s --output %s"%(args.cuts,args.outdir)
print comm
os.system(comm)

print "Priors for this sample:"
print dz
print frac

#Ddz = np.array([0.0153, 0.0130, 0.0112, 0.0136]) / np.sqrt(frac) 
Ddz = np.loadtxt("%s/cosmos-prior_widths.txt"%args.outdir)

export_priors(dz,Ddz, "%s/mcal_priors.ini"%args.outdir)
