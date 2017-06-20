import numpy as np
import pylab as plt
import fitsio, yaml, os
import argparse
import socket
plt.switch_backend("pdf")
plt.style.use("y1a1")

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--cuts', default='none', action='store')
parser.add_argument('--output', default="plots", action='store')
args = parser.parse_args()

hostname = socket.gethostname()
if 'fornax' in hostname:
    path="/share/des/disc4/samuroff/ia/pz_calibration/"
else:
    path="/global/cscratch1/sd/sws/y1/photoz/"
f="%s/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.fits"%path


def wavg(x,wt): return sum(x*wt)/sum(wt)


def get_mask(lower,upper,col_name,catalog):
    return (catalog[col_name]<upper) & (catalog[col_name]>lower)


# (5) get uncertainty due to morphology matching
# we have two versions of the COSMOS resampled/reweighted catalogs:
# one where we only try to match fluxes (=colors/magnitudes), v3, and one where we also try to match 
# the intrinsic sizes of galaxies, v4
# The latter is our fiducial for calibration. Lensing selection does depend on galaxy size quite a bit.
# But (i) our size matching is probably not perfect 
# (in fact, we don't even try to match sizes for ~30 per-cent of galaxies because it would make the matched colors too much worse)
# and (ii) size is only one proxy (that we measure and therefore can match) of whether a galaxy makes it into the shape catalog, and with what weight
# higher-order effects might matter as well, at a lower level
# We are conservative and back-of-the-envelope on this: half the mean z difference of the v3 and v4 sample 
# is assigned as a systematic uncertainty due to galaxy morphology matching

# v4 catalogs (with size matching)
f="%s/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.fits"%path
rphot=fitsio.read(f) # resampled photometry, 200.000 random objects from the science sample
rbpz=fitsio.read("%s/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.BPZ.fits"%path) # BPZ run on that MOF photometry

# v3 catalogs (without size matching)
rphot3=fitsio.read("%s/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v3.fits"%path)
rbpz3=fitsio.read("%s/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v3.BPZ.fits"%path)


if args.cuts!='none':
    cuts = yaml.load(open(args.cuts))
else:
   cuts = {}
mcal_mask = np.ones(rbpz.size).astype(bool)
mcal_mask3 = np.ones(rbpz3.size).astype(bool)
lower = -1e5
upper = 1e5
col = "TEMPLATE_TYPE"
for col in cuts.keys():
    print "Cutting %s"%col
    lower,upper,dt = cuts[col].split()
    exec("lower=%s(lower)"%dt)
    exec("upper=%s(upper)"%dt)
    mcal_mask = mcal_mask & (rbpz[col]>lower) & (rbpz[col]<upper)
    mcal_mask3 = mcal_mask3 & (rbpz3[col]>lower) & (rbpz3[col]<upper)


# assigning to redshift bins, v4
lmask0,=np.where((rbpz['MEAN_Z']>0.20)&(rbpz['MEAN_Z']<1.30) & mcal_mask)
lmask1,=np.where((rbpz['MEAN_Z']>0.20)&(rbpz['MEAN_Z']<0.43) & mcal_mask)
lmask2,=np.where((rbpz['MEAN_Z']>0.43)&(rbpz['MEAN_Z']<0.63) & mcal_mask)
lmask3,=np.where((rbpz['MEAN_Z']>0.63)&(rbpz['MEAN_Z']<0.90) & mcal_mask)
lmask4,=np.where((rbpz['MEAN_Z']>0.90)&(rbpz['MEAN_Z']<1.30) & mcal_mask)

# assigning to redshift bins, v3
l3mask0,=np.where((rbpz3['MEAN_Z']>0.20)&(rbpz3['MEAN_Z']<1.30) & mcal_mask3)
l3mask1,=np.where((rbpz3['MEAN_Z']>0.20)&(rbpz3['MEAN_Z']<0.43) & mcal_mask3)
l3mask2,=np.where((rbpz3['MEAN_Z']>0.43)&(rbpz3['MEAN_Z']<0.63) & mcal_mask3)
l3mask3,=np.where((rbpz3['MEAN_Z']>0.63)&(rbpz3['MEAN_Z']<0.90) & mcal_mask3)
l3mask4,=np.where((rbpz3['MEAN_Z']>0.90)&(rbpz3['MEAN_Z']<1.30) & mcal_mask3)

# mean z differences
diff=[]
print("bin 0 difference=",(wavg(rbpz['redshift'][lmask0],rphot['R11'][lmask0]+rphot['R22'][lmask0])
      -wavg(rbpz3['redshift'][l3mask0],rphot3['R11'][l3mask0]+rphot3['R22'][l3mask0]))/2)

print("bin 1 difference=",(wavg(rbpz['redshift'][lmask1],rphot['R11'][lmask1]+rphot['R22'][lmask1])
      -wavg(rbpz3['redshift'][l3mask1],rphot3['R11'][l3mask1]+rphot3['R22'][l3mask1]))/2)
diff.append((wavg(rbpz['redshift'][lmask1],rphot['R11'][lmask1]+rphot['R22'][lmask1])
      -wavg(rbpz3['redshift'][l3mask1],rphot3['R11'][l3mask1]+rphot3['R22'][l3mask1]))/2)
print("bin 2 difference=",(wavg(rbpz['redshift'][lmask2],rphot['R11'][lmask2]+rphot['R22'][lmask2])
      -wavg(rbpz3['redshift'][l3mask2],rphot3['R11'][l3mask2]+rphot3['R22'][l3mask2]))/2)
diff.append((wavg(rbpz['redshift'][lmask2],rphot['R11'][lmask2]+rphot['R22'][lmask2])
      -wavg(rbpz3['redshift'][l3mask2],rphot3['R11'][l3mask2]+rphot3['R22'][l3mask2]))/2)
print("bin 3 difference=",(wavg(rbpz['redshift'][lmask3],rphot['R11'][lmask3]+rphot['R22'][lmask3])
      -wavg(rbpz3['redshift'][l3mask3],rphot3['R11'][l3mask3]+rphot3['R22'][l3mask3]))/2)
diff.append((wavg(rbpz['redshift'][lmask3],rphot['R11'][lmask3]+rphot['R22'][lmask3])
      -wavg(rbpz3['redshift'][l3mask3],rphot3['R11'][l3mask3]+rphot3['R22'][l3mask3]))/2)
print("bin 4 difference=",(wavg(rbpz['redshift'][lmask4],rphot['R11'][lmask4]+rphot['R22'][lmask4])
      -wavg(rbpz3['redshift'][l3mask4],rphot3['R11'][l3mask4]+rphot3['R22'][l3mask4]))/2)
diff.append((wavg(rbpz['redshift'][lmask4],rphot['R11'][lmask4]+rphot['R22'][lmask4])
      -wavg(rbpz3['redshift'][l3mask4],rphot3['R11'][l3mask4]+rphot3['R22'][l3mask4]))/2)

diff=np.array(diff)
print("rms difference=",np.sqrt(sum(diff*diff)/len(diff)))
dz_morphology = np.array([np.sqrt(sum(diff*diff)/len(diff))] * 4)
print("The rms difference goes to https://docs.google.com/document/d/1Bo_zMI1S2F-Han7KkxAS-tDHErKJq0nt9YPu9Vl6EKc")
print("\"COSMOS resampling morphology uncertainty\"")



# Try calculating a covariance matrix
ngal = rbpz['MEAN_Z'].size
n0=ngal/5
indices = np.random.choice(ngal,ngal, replace=False)
vec=[[],[],[],[]]
for i in xrange(5):
    subset = indices[(i*n0):(i+1)*n0]
    lmask0,=np.where((rbpz['MEAN_Z']>0.20)&(rbpz['MEAN_Z']<1.30) & mcal_mask)
    lmask1,=np.where((rbpz['MEAN_Z']>0.20)&(rbpz['MEAN_Z']<0.43) & mcal_mask)
    lmask2,=np.where((rbpz['MEAN_Z']>0.43)&(rbpz['MEAN_Z']<0.63) & mcal_mask)
    lmask3,=np.where((rbpz['MEAN_Z']>0.63)&(rbpz['MEAN_Z']<0.90) & mcal_mask)
    lmask4,=np.where((rbpz['MEAN_Z']>0.90)&(rbpz['MEAN_Z']<1.30) & mcal_mask)
    l3mask0,=np.where((rbpz3['MEAN_Z']>0.20)&(rbpz3['MEAN_Z']<1.30) & mcal_mask3)
    l3mask1,=np.where((rbpz3['MEAN_Z']>0.20)&(rbpz3['MEAN_Z']<0.43) & mcal_mask3)
    l3mask2,=np.where((rbpz3['MEAN_Z']>0.43)&(rbpz3['MEAN_Z']<0.63) & mcal_mask3)
    l3mask3,=np.where((rbpz3['MEAN_Z']>0.63)&(rbpz3['MEAN_Z']<0.90) & mcal_mask3)
    l3mask4,=np.where((rbpz3['MEAN_Z']>0.90)&(rbpz3['MEAN_Z']<1.30) & mcal_mask3)

    l0 = lmask0[np.in1d(lmask0,subset)]
    l1 = lmask1[np.in1d(lmask1,subset)]
    l2 = lmask2[np.in1d(lmask2,subset)]
    l3 = lmask3[np.in1d(lmask3,subset)]
    l4 = lmask4[np.in1d(lmask4,subset)]
    l30 = l3mask0[np.in1d(l3mask0,subset)]
    l31 = l3mask1[np.in1d(l3mask1,subset)]
    l32 = l3mask2[np.in1d(l3mask2,subset)]
    l33 = l3mask3[np.in1d(l3mask3,subset)]
    l34 = l3mask4[np.in1d(l3mask4,subset)]
    #import pdb ; pdb.set_trace()
    vec[0].append((wavg(rbpz['redshift'][l1],rphot['R11'][l1]+rphot['R22'][l1])-wavg(rbpz3['redshift'][l31],rphot3['R11'][l31]+rphot3['R22'][l31]))/2)
    vec[1].append((wavg(rbpz['redshift'][l2],rphot['R11'][l2]+rphot['R22'][l2])-wavg(rbpz3['redshift'][l32],rphot3['R11'][l32]+rphot3['R22'][l32]))/2)
    vec[2].append((wavg(rbpz['redshift'][l3],rphot['R11'][l3]+rphot['R22'][l3])-wavg(rbpz3['redshift'][l33],rphot3['R11'][l33]+rphot3['R22'][l33]))/2)
    vec[3].append((wavg(rbpz['redshift'][l4],rphot['R11'][l4]+rphot['R22'][l4])-wavg(rbpz3['redshift'][l34],rphot3['R11'][l34]+rphot3['R22'][l34]))/2)
    

cov=np.cov(vec)
import pdb ; pdb.set_trace()
