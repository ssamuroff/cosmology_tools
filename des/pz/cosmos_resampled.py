import numpy as np
import pylab as plt
import fitsio, yaml, os
import argparse
plt.switch_backend("pdf")
plt.style.use("y1a1")

print "----------------------------------------"
print "n(z) recalibration and diagnostic script"
print "Courtesy of Daniel Gruen"
print "-----------------------------------------"

print "I've hardcoded some of the paths for now, but this is fine on fornax."


parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--cuts', default='none', action='store')
parser.add_argument('--output', default="plots", action='store')
parser.add_argument('--width', action='store_true')
args = parser.parse_args()

# (2) here begins the real deal: use these files to get photo-z metrics
# we have a single resampling with color+magnitude+size matching and without a photometric offset applied
# this is in cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.fits

f="/share/des/disc4/samuroff/ia/pz_calibration/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.fits"
rphot=fitsio.read(f) # resampled photometry, 200.000 random objects from the science sample
rmofbpz=fitsio.read("/share/des/disc4/samuroff/ia/pz_calibration/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.BPZMOF.fits") # BPZ run on that MOF photometry
rmcalbpz=fitsio.read("/share/des/disc4/samuroff/ia/pz_calibration/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.BPZ.fits") # BPZ run on that MCAL photometry

print(
np.max(rphot['coadd_objects_id']-rmofbpz['coadd_objects_id']),
np.max(rphot['coadd_objects_id']-rmcalbpz['coadd_objects_id'])) # check IDs are the same


if args.cuts!='none':
    cuts = yaml.load(open(args.cuts))
else:
   cuts = {}
mcal_mask = np.ones(rmcalbpz.size).astype(bool)
mof_mask = np.ones(rmofbpz.size).astype(bool)
for col in cuts.keys():
    print "Cutting %s"%col
    lower,upper,dt = cuts[col].split()
    exec("lower=%s(lower)"%dt)
    exec("upper=%s(upper)"%dt)
    mcal_mask = mcal_mask & (rmcalbpz[col]>lower) & (rmcalbpz[col]<upper)
    mof_mask = mof_mask & (rmofbpz[col]>lower) & (rmofbpz[col]<upper)

os.system("mkdir -p %s"%args.output)

# (2a) look at the photometric chi^2 distribution of the matches,
# i.e., how different are the matched COSMOS galaxies to the underlying science galaxies, 
#based on their flux measurements and errors
plt.hist(rphot['matched_chi2'][mof_mask],bins=100,range=(0,20))
plt.plot((4,4),(0,10000),color='black') # 4 degrees of freedom
plt.ylim(0,7000)
plt.xlabel(r"$\chi^2$")
#plt.savefig("%s/chi2-cosmos_matches.pdf"%args.output)
plt.close()

# (2b) our resampling algorithm only works for COSMOS galaxies with smaller flux errors
# so check how many COSMOS galaxies have had smaller errors than the metacal galaxy in each case
x=plt.hist(rphot["nerrmask"][mof_mask],bins=100)
plt.xlim(0,1.35e5)
plt.xlabel(r"number of COSMOS galaxies w/ $<$ errors")
#plt.savefig("%s/num-cosmos-low_flux_error.pdf"%args.output)
plt.close()

print(sum(rphot['nsizemask'][mof_mask]>1)*1.0/len(rphot[mof_mask])*100,"percent of galaxies are matched in size as well as fluxes")
print("(for the rest, matching in size makes photometric chi^2 worse by more than 4)")



# (3) get mean true redshift of each MCAL BPZ mean_z bin from COSMOS

# assign to bins in MCAL MEAN_Z
# 0=all: 0.2<mean_z<1.3
# 1=0.2..0.43
# 2=0.43..0.63
# 3=0.63..0.90
# 4=0.90..1.30
lmask0,=np.where((rmcalbpz['MEAN_Z']>0.20)&(rmcalbpz['MEAN_Z']<1.30) & mcal_mask)
lmask1,=np.where((rmcalbpz['MEAN_Z']>0.20)&(rmcalbpz['MEAN_Z']<0.43) & mcal_mask)
lmask2,=np.where((rmcalbpz['MEAN_Z']>0.43)&(rmcalbpz['MEAN_Z']<0.63) & mcal_mask)
lmask3,=np.where((rmcalbpz['MEAN_Z']>0.63)&(rmcalbpz['MEAN_Z']<0.90) & mcal_mask)
lmask4,=np.where((rmcalbpz['MEAN_Z']>0.90)&(rmcalbpz['MEAN_Z']<1.30) & mcal_mask)

def wavg(x,wt):
    return sum(x*wt)/sum(wt)

print("all <z>=",wavg(rphot['redshift'][lmask0],rphot['R11'][lmask0]+rphot['R22'][lmask0]))
print("bin 1 <z>=",wavg(rphot['redshift'][lmask1],rphot['R11'][lmask1]+rphot['R22'][lmask1]))
print("bin 2 <z>=",wavg(rphot['redshift'][lmask2],rphot['R11'][lmask2]+rphot['R22'][lmask2]))
print("bin 3 <z>=",wavg(rphot['redshift'][lmask3],rphot['R11'][lmask3]+rphot['R22'][lmask3]))
print("bin 4 <z>=",wavg(rphot['redshift'][lmask4],rphot['R11'][lmask4]+rphot['R22'][lmask4]))
print("these numbers go to https://docs.google.com/document/d/1Bo_zMI1S2F-Han7KkxAS-tDHErKJq0nt9YPu9Vl6EKc")
print("\"COSMOS <z>, metacal\"")

#zavg = [wavg(rphot['redshift'][lm],rphot['R11'][lm]+rphot['R22'][lm]) for lm in [lmask0,lmask1,lmask2,lmask3,lmask4]]
#np.savetxt("%s/cosmos-zmean-bins.txt"%args.output,zavg)

# get Troxel his histograms
lmask0p = (rphot['R11'][lmask0]+rphot['R22'][lmask0]>0)
lmask0m = (rphot['R11'][lmask0]+rphot['R22'][lmask0]<=0)
lmask1p = (rphot['R11'][lmask1]+rphot['R22'][lmask1]>0)
lmask1m = (rphot['R11'][lmask1]+rphot['R22'][lmask1]<=0)
lmask2p = (rphot['R11'][lmask2]+rphot['R22'][lmask2]>0)
lmask2m = (rphot['R11'][lmask2]+rphot['R22'][lmask2]<=0)
lmask3p = (rphot['R11'][lmask3]+rphot['R22'][lmask3]>0)
lmask3m = (rphot['R11'][lmask3]+rphot['R22'][lmask3]<=0)
lmask4p = (rphot['R11'][lmask4]+rphot['R22'][lmask4]>0)
lmask4m = (rphot['R11'][lmask4]+rphot['R22'][lmask4]<=0)


p0=plt.hist(rphot['redshift'][lmask0][lmask0p],weights=rphot['R11'][lmask0][lmask0p]+rphot['R22'][lmask0][lmask0p],bins=400,range=(0,4))
p1=plt.hist(rphot['redshift'][lmask1][lmask1p],weights=rphot['R11'][lmask1][lmask1p]+rphot['R22'][lmask1][lmask1p],bins=400,range=(0,4),alpha=0.5)
p2=plt.hist(rphot['redshift'][lmask2][lmask2p],weights=rphot['R11'][lmask2][lmask2p]+rphot['R22'][lmask2][lmask2p],bins=400,range=(0,4),alpha=0.5)
p3=plt.hist(rphot['redshift'][lmask3][lmask3p],weights=rphot['R11'][lmask3][lmask3p]+rphot['R22'][lmask3][lmask3p],bins=400,range=(0,4),alpha=0.5)
p4=plt.hist(rphot['redshift'][lmask4][lmask4p],weights=rphot['R11'][lmask4][lmask4p]+rphot['R22'][lmask4][lmask4p],bins=400,range=(0,4),alpha=0.5)
plt.figure()
m0=plt.hist(rphot['redshift'][lmask0][lmask0m],weights=-rphot['R11'][lmask0][lmask0m]-rphot['R22'][lmask0][lmask0m],bins=400,range=(0,4))
m1=plt.hist(rphot['redshift'][lmask1][lmask1m],weights=-rphot['R11'][lmask1][lmask1m]-rphot['R22'][lmask1][lmask1m],bins=400,range=(0,4),alpha=0.5)
m2=plt.hist(rphot['redshift'][lmask2][lmask2m],weights=-rphot['R11'][lmask2][lmask2m]-rphot['R22'][lmask2][lmask2m],bins=400,range=(0,4),alpha=0.5)
m3=plt.hist(rphot['redshift'][lmask3][lmask3m],weights=-rphot['R11'][lmask3][lmask3m]-rphot['R22'][lmask3][lmask3m],bins=400,range=(0,4),alpha=0.5)
m4=plt.hist(rphot['redshift'][lmask4][lmask4m],weights=-rphot['R11'][lmask4][lmask4m]-rphot['R22'][lmask4][lmask4m],bins=400,range=(0,4),alpha=0.5)

np.savetxt("%s/bin0_metacal_COSMOS.tab"%args.output,np.array([p0[1][:-1],p0[1][1:],p0[0]+m0[0]]).transpose(),header="zmin zmax weight")
np.savetxt("%s/bin1_metacal_COSMOS.tab"%args.output,np.array([p1[1][:-1],p1[1][1:],p1[0]+m1[0]]).transpose(),header="zmin zmax weight")
np.savetxt("%s/bin2_metacal_COSMOS.tab"%args.output,np.array([p2[1][:-1],p2[1][1:],p2[0]+m2[0]]).transpose(),header="zmin zmax weight")
np.savetxt("%s/bin3_metacal_COSMOS.tab"%args.output,np.array([p3[1][:-1],p3[1][1:],p3[0]+m3[0]]).transpose(),header="zmin zmax weight")
np.savetxt("%s/bin4_metacal_COSMOS.tab"%args.output,np.array([p4[1][:-1],p4[1][1:],p4[0]+m4[0]]).transpose(),header="zmin zmax weight")


# (3) get mean true redshift of each MOF BPZ mean_z bin, im3shape weighted, from COSMOS



#f="wl_class.im3shape.METACAL_MOF.cosmos.v3.fits"
#rphot=fitsio.read(f) # resampled photometry, 200.000 random objects from the science sample
#rmofbpz=fitsio.read("/share/des/disc4/samuroff/ia/pz_calibration/cats/wl_class.im3shape.METACAL_MOF.cosmos.v3.BPZMOF.fits") # BPZ run on that MOF photometry, now with MOF prior

# assign to bins in MOF MEAN_Z
# 0=all: 0.2<mean_z<1.3
# 1=0.2..0.43
# 2=0.43..0.63
# 3=0.63..0.90
# 4=0.90..1.30
#lmask0,=np.where((rmofbpz['MEAN_Z']>0.20)&(rmofbpz['MEAN_Z']<1.30))
#lmask1,=np.where((rmofbpz['MEAN_Z']>0.20)&(rmofbpz['MEAN_Z']<0.43))
#lmask2,=np.where((rmofbpz['MEAN_Z']>0.43)&(rmofbpz['MEAN_Z']<0.63))
#lmask3,=np.where((rmofbpz['MEAN_Z']>0.63)&(rmofbpz['MEAN_Z']<0.90))
#lmask4,=np.where((rmofbpz['MEAN_Z']>0.90)&(rmofbpz['MEAN_Z']<1.30))

#print("all <z>=",wavg(rphot['redshift'][lmask0],rphot['weight'][lmask0]*(1.+rphot['m'][lmask0])))
#print("bin 1 <z>=",wavg(rphot['redshift'][lmask1],rphot['weight'][lmask1]*(1.+rphot['m'][lmask1])))
#print("bin 2 <z>=",wavg(rphot['redshift'][lmask2],rphot['weight'][lmask2]*(1.+rphot['m'][lmask2])))
#print("bin 3 <z>=",wavg(rphot['redshift'][lmask3],rphot['weight'][lmask3]*(1.+rphot['m'][lmask3])))
#print("bin 4 <z>=",wavg(rphot['redshift'][lmask4],rphot['weight'][lmask4]*(1.+rphot['m'][lmask4])))
#print("these numbers go to https://docs.google.com/document/d/1Bo_zMI1S2F-Han7KkxAS-tDHErKJq0nt9YPu9Vl6EKc")
#print("\"COSMOS <z>, im3shape\"")

# (repeat 3) get mean true redshift of each MCAL BPZ mean_z bin from COSMOS, 
# this time with a BPZ prior based on MOF mag i

f="/share/des/disc4/samuroff/ia/pz_calibration/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.fits"
rphot=fitsio.read(f) # resampled photometry, 200.000 random objects from the science sample
rmcalbpz=fitsio.read("/share/des/disc4/samuroff/ia/pz_calibration/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.BPZ.fits") # BPZ run on that MCAL photometry

mcal_mask = np.ones(rmcalbpz.size).astype(bool)
mof_mask = np.ones(rmofbpz.size).astype(bool)
for col in cuts.keys():
    print "Cutting %s"%col
    lower,upper,dt = cuts[col].split()
    exec("lower=%s(lower)"%dt)
    exec("upper=%s(upper)"%dt)
    mcal_mask = mcal_mask & (rmcalbpz[col]>lower) & (rmcalbpz[col]<upper)
    mof_mask = mof_mask & (rmofbpz[col]>lower) & (rmofbpz[col]<upper)

# assign to bins in MCAL MEAN_Z
# 0=all: 0.2<mean_z<1.3
# 1=0.2..0.43
# 2=0.43..0.63
# 3=0.63..0.90
# 4=0.90..1.30
lmask0,=np.where((rmcalbpz['MEAN_Z']>0.20)&(rmcalbpz['MEAN_Z']<1.30) & mcal_mask)
lmask1,=np.where((rmcalbpz['MEAN_Z']>0.20)&(rmcalbpz['MEAN_Z']<0.43) & mcal_mask)
lmask2,=np.where((rmcalbpz['MEAN_Z']>0.43)&(rmcalbpz['MEAN_Z']<0.63) & mcal_mask)
lmask3,=np.where((rmcalbpz['MEAN_Z']>0.63)&(rmcalbpz['MEAN_Z']<0.90) & mcal_mask)
lmask4,=np.where((rmcalbpz['MEAN_Z']>0.90)&(rmcalbpz['MEAN_Z']<1.30) & mcal_mask)

print("what we did before first: MCAL mag i for prior")
print("all <z>=",wavg(rphot['redshift'][lmask0],rphot['R11'][lmask0]+rphot['R22'][lmask0]))
print("bin 1 <z>=",wavg(rphot['redshift'][lmask1],rphot['R11'][lmask1]+rphot['R22'][lmask1]))
print("bin 2 <z>=",wavg(rphot['redshift'][lmask2],rphot['R11'][lmask2]+rphot['R22'][lmask2]))
print("bin 3 <z>=",wavg(rphot['redshift'][lmask3],rphot['R11'][lmask3]+rphot['R22'][lmask3]))
print("bin 4 <z>=",wavg(rphot['redshift'][lmask4],rphot['R11'][lmask4]+rphot['R22'][lmask4]))

zavg = [wavg(rphot['redshift'][lm],rphot['R11'][lm]+rphot['R22'][lm]) for lm in [lmask0,lmask1,lmask2,lmask3,lmask4]]
np.savetxt("%s/cosmos-zmean-bins-mcal_iprior.txt"%args.output,zavg)

f="/share/des/disc4/samuroff/ia/pz_calibration/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.fits"
rphot=fitsio.read(f) # resampled photometry, 200.000 random objects from the science sample
rmcalbpz=fitsio.read("/share/des/disc4/samuroff/ia/pz_calibration/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.BPZ.MCALWITHMOFPRIOR.fits") 
# sorry, that catalog is not part of the gz file
# you can get it straight from 
# https://www.slac.stanford.edu/~dgruen/y1cat/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.BPZ.MCALWITHMOFPRIOR.fits
# BPZ run on that MCAL photometry, this time using MOF i for the prior

# assign to bins in MCAL MEAN_Z with MOF i prior
# 0=all: 0.2<mean_z<1.3
# 1=0.2..0.43
# 2=0.43..0.63
# 3=0.63..0.90
# 4=0.90..1.30

mcal_mask = np.ones(rmcalbpz.size).astype(bool)
mof_mask = np.ones(rmofbpz.size).astype(bool)
for col in cuts.keys():
    print "Cutting %s"%col
    lower,upper,dt = cuts[col].split()
    exec("lower=%s(lower)"%dt)
    exec("upper=%s(upper)"%dt)
    mcal_mask = mcal_mask & (rmcalbpz[col]>lower) & (rmcalbpz[col]<upper)
    mof_mask = mof_mask & (rmofbpz[col]>lower) & (rmofbpz[col]<upper)


lmask0,=np.where((rmcalbpz['MEAN_Z']>0.20)&(rmcalbpz['MEAN_Z']<1.30) & mcal_mask)
lmask1,=np.where((rmcalbpz['MEAN_Z']>0.20)&(rmcalbpz['MEAN_Z']<0.43) & mcal_mask)
lmask2,=np.where((rmcalbpz['MEAN_Z']>0.43)&(rmcalbpz['MEAN_Z']<0.63) & mcal_mask)
lmask3,=np.where((rmcalbpz['MEAN_Z']>0.63)&(rmcalbpz['MEAN_Z']<0.90) & mcal_mask)
lmask4,=np.where((rmcalbpz['MEAN_Z']>0.90)&(rmcalbpz['MEAN_Z']<1.30) & mcal_mask)

print("what Ben did: MOF mag i for prior")
print("all <z>=",wavg(rphot['redshift'][lmask0],rphot['R11'][lmask0]+rphot['R22'][lmask0]))
print("bin 1 <z>=",wavg(rphot['redshift'][lmask1],rphot['R11'][lmask1]+rphot['R22'][lmask1]))
print("bin 2 <z>=",wavg(rphot['redshift'][lmask2],rphot['R11'][lmask2]+rphot['R22'][lmask2]))
print("bin 3 <z>=",wavg(rphot['redshift'][lmask3],rphot['R11'][lmask3]+rphot['R22'][lmask3]))
print("bin 4 <z>=",wavg(rphot['redshift'][lmask4],rphot['R11'][lmask4]+rphot['R22'][lmask4]))
print("these numbers go to https://docs.google.com/document/d/1Bo_zMI1S2F-Han7KkxAS-tDHErKJq0nt9YPu9Vl6EKc")
print("\"COSMOS <z>, metacal\"")

zavg = [wavg(rphot['redshift'][lm],rphot['R11'][lm]+rphot['R22'][lm]) for lm in [lmask0,lmask1,lmask2,lmask3,lmask4]]
np.savetxt("%s/cosmos-zmean-bins-mof_iprior.txt"%args.output,zavg)

if not args.width: exit()
