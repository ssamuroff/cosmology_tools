import numpy as np
import pylab as plt
import fitsio, yaml, os
import argparse
plt.switch_backend("pdf")
plt.style.use("y1a1")


print "Calculating prior widths --  this may take a few minutes to finish."

# There are five error components here
# Calculate them and combine them in quadrature

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--cuts', default='none', action='store')
parser.add_argument('--output', default="plots", action='store')
args = parser.parse_args()

# (4) get uncertainty due to size of the resampled catalog from a bootstrap estimation
path = "/share/des/disc4/samuroff/ia/pz_calibration/"
f="%s/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.fits"%path
rphot=fitsio.read(f) # resampled photometry, 200.000 random objects from the science sample
rmcalbpz=fitsio.read("%s/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.BPZ.fits"%path) # BPZ run on that MCAL photometry
rmofbpz=fitsio.read("%s/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.BPZMOF.fits"%path) # BPZ run on that MOF photometry

if args.cuts!='none':
    cuts = yaml.load(open(args.cuts))
else:
   cuts = {}
mcal_mask = np.ones(rmcalbpz.size).astype(bool)
mof_mask = np.ones(rmofbpz.size).astype(bool)
lower = -1e5
upper = 1e5
col = "TEMPLATE_TYPE"
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


# redshift bins based on METACAL BPZ mean-z estimate
l0=[] # all: 0.2<mean_z<1.3
l1=[] # 0.2..0.43
l2=[] # 0.43..0.63
l3=[] # 0.63..0.90
l4=[] # 0.90..1.30


def wavg(x,wt): return sum(x*wt)/sum(wt)


def get_mask(lower,upper,col_name,catalog):
    return (catalog[col_name]<upper) & (catalog[col_name]>lower)

smask = np.ones(rmcalbpz.size).astype(bool)
for i in range(500): # do 500 bootstrap resamplings; 
# this is going to take a while, but if you get impatient you can interrupt the kernel and still go on 
    np.random.seed(i) # different random seed each time
    sphot=np.random.choice(rphot, size=len(rphot)) # pick a random sample with repetition
    np.random.seed(i) # same random seed as above
    sbpz=np.random.choice(rmcalbpz, size=len(rphot)) # pick same random sample from METACAL BPZ catalog for binning
    smask = get_mask(lower,upper,col,sbpz)
    
    # assign to the 4+1 bins
    lmask0,=np.where((sbpz['MEAN_Z']>0.20)&(sbpz['MEAN_Z']<1.30) & smask)
    lmask1,=np.where((sbpz['MEAN_Z']>0.20)&(sbpz['MEAN_Z']<0.43) & smask)
    lmask2,=np.where((sbpz['MEAN_Z']>0.43)&(sbpz['MEAN_Z']<0.63) & smask)
    lmask3,=np.where((sbpz['MEAN_Z']>0.63)&(sbpz['MEAN_Z']<0.90) & smask)
    lmask4,=np.where((sbpz['MEAN_Z']>0.90)&(sbpz['MEAN_Z']<1.30) & smask)
    
    for l,m in zip([l0,l1,l2,l3,l4],[lmask0,lmask1,lmask2,lmask3,lmask4]):
        # save mean true-z of each bin
        l.append(wavg(sbpz['redshift'][m],sphot['R11'][m]+sphot['R22'][m]))


# convert lists to numpy arrays so we can better handle them with numpy functions
l0=np.array(l0)
l1=np.array(l1)
l2=np.array(l2)
l3=np.array(l3)
l4=np.array(l4)

dz_stat = []

for La,a in zip([l0,l1,l2,l3,l4],range(5)):
    for Lb,b in zip([l0,l1,l2,l3,l4],range(5)):
        print("COV(z",a,", z",b,")=",np.sum(La*Lb)/len(La)-np.sum(La)*sum(Lb)/len(La)/len(Lb))
        if(a==b):
            print("stdev z",a,"=",np.sqrt(np.sum(La*Lb)/len(La)-np.sum(La)*sum(Lb)/len(La)/len(Lb)))
            dz_stat.append(np.sqrt(np.sum(La*Lb)/len(La)-np.sum(La)*sum(Lb)/len(La)/len(Lb)))

print("the lines with stdev go into https://docs.google.com/document/d/1Bo_zMI1S2F-Han7KkxAS-tDHErKJq0nt9YPu9Vl6EKc")
print("as \"COSMOS resampling stat. uncertainty\"")
print("the fact that the off-diagonal components of the COV are much smaller than the diagonal ones tells us that these are uncorrelated between bins (as they should be when bins don't overlap)")



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




# (6) have a look at COSMOS resamplings in the Buzzard simulations

# to this end, we have used a lensing-sample-like sample of galaxies in the Buzzard simulations
# we cut out COSMOS-shaped footprints from the Buzzard simulations and find a best 
# match to 200.000 galaxies from the lensing sample the same way as we do in the data
# notable difference: there is no relative weighting of the shape sources

b = [1, 2, 3, 5, 6, 8, 9, 10, 11, 13, 14, 15, 16, 17, 18, 20, 21, 22, 24, 25, 26, 27, 28, 29, 30, 31, 33, 34, 35, 36, 37, 38, 40, 41, 42, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 67, 68, 69, 70, 71, 72, 73, 75, 76, 77, 78, 79, 80, 81, 82, 85, 88, 89, 90, 91, 93, 97, 98, 100, 101, 103, 105, 107, 108, 109, 111, 112, 113, 114, 115, 117, 118, 119, 121, 122, 123, 124, 125, 126, 128, 129, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 149, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 167, 169, 171, 172, 173, 174, 176, 177, 178, 180, 181, 182, 184, 185, 187, 188, 190, 192, 194, 195, 196, 197, 198, 199, 200, 201, 202, 204, 205, 206, 208, 209, 212, 213, 214, 215, 216, 217, 218, 219, 222, 224, 226, 227, 228, 229, 230, 231, 232, 233, 235, 236, 238, 239, 240, 241, 243, 244, 245, 246, 247, 248, 249, 250, 252, 253, 254, 255, 256, 257, 258, 261, 262, 263, 264, 265, 266, 267, 268, 269, 271, 272, 274, 275, 277, 278, 279, 280, 281, 283, 284, 286, 287, 288, 289, 290, 291, 293, 294, 295, 298, 299, 301, 303, 305, 307, 308, 309, 310, 312, 313, 316, 317, 318, 319, 320, 321, 322, 323, 324, 325, 329, 330, 331, 332, 334, 335, 336, 337, 339, 340, 342, 343, 344, 349, 350, 351, 353, 356, 357, 358, 359, 361, 362, 363, 364, 365, 366, 368, 372, 373, 374, 376, 378, 379, 380, 381, 382, 383, 384, 385, 386, 387, 388, 389, 390, 391, 392, 394, 395, 397, 399, 400, 401, 402, 404, 405, 406, 408, 409, 410, 411, 412, 413, 414, 416, 417, 419, 420, 421, 425, 426, 429, 432, 434, 435, 437, 441, 442, 444, 447, 448, 449, 450, 452, 453, 457, 458, 460, 461, 462, 463, 466, 467, 469, 471, 473, 474, 476, 477, 479, 480, 481, 483, 484, 485, 486, 487, 488, 491, 492, 493, 494, 495, 496, 497, 498]
print("using",len(b),"good COSMOS resamplings")
# these are the IDs of good COSMOS footprints in Buzzard

d=fitsio.read("%s/cats/buzzard/cosmos_resampled_1.fits"%path)


# ok, let's get to cosmic variance (7) and systematic uncertainty (8)
# we get COSMIC variance from the scatter of mean true z in the resamplings of the 368 simulated COSMOS footprints
# we get some idea of systematic uncertainty just from the difference of true redshift of the lensing sample 
#    and the resampled-matched COSMOS sample
# for this, everything is binned by BPZ mean-z estimated from the resampled photometry, as always


# organized as Python dictionaries because it's convenient here
bpzmeans  = {}   # means of BPZ estimate of resampled COSMOS
truemeans_c = {} # means of true z of resampled COSMOS
truemeans_l = {} # means of true z of underlying lensing sample that the COSMOS catalog was matched to

for i in range(5):
    bpzmeans[i] = []
    truemeans_c[i] = []
    truemeans_l[i] = []

dmask = np.ones(rmcalbpz.size).astype(bool)

for c in b:
    d = fitsio.read("%s/cats/buzzard/cosmos_resampled_"%(path) + str(c)+".fits")
    dbpz = fitsio.read("%s/cats/buzzard/cosmos_resampled_"%(path) + str(c)+".BPZ.fits") 
    dmask = dmask & get_mask(lower,upper,col,dbpz)
    # BPZ run on the resampled photometry, for binning
    
    assert(np.max(np.abs(dbpz['coadd_objects_id']-d['coadd_objects_id']))==0) # assert same ordering
    
    # make bins, as you know from earlier
    bins={}
    bins[0] ,= np.where((dbpz['MEAN_Z']>0.2)&(dbpz['MEAN_Z']<=1.3) & dmask)
    bins[1] ,= np.where((dbpz['MEAN_Z']>0.2)&(dbpz['MEAN_Z']<=0.43) & dmask)
    bins[2] ,= np.where((dbpz['MEAN_Z']>0.43)&(dbpz['MEAN_Z']<=0.63) & dmask)
    bins[3] ,= np.where((dbpz['MEAN_Z']>0.63)&(dbpz['MEAN_Z']<=0.90) & dmask)
    bins[4] ,= np.where((dbpz['MEAN_Z']>0.90)&(dbpz['MEAN_Z']<=1.3) & dmask)

    for i in range(5):
        bpzmeans[i].append(np.average(dbpz['Z_MC'][bins[i]]))
        truemeans_c[i].append(np.average(d['redshift'][bins[i]]))
        truemeans_l[i].append(np.average(d['buzzard_redshift'][bins[i]]))
 
dz_sys = np.sqrt(np.mean(np.array([np.mean(np.array(truemeans_l[i])-np.array(truemeans_c[i])) for i in [1,2,3,4]])**2))
dz_sys = [dz_sys] * 4
dz_cv = np.array([np.std(truemeans_c[i]) for i in [1,2,3,4]])

# finally, let's do uncertainty from photometric calibration (9) 
# we get this from the scatter of mean <Z_MC> among 200 resampling/reweightings of COSMOS, 
# where we have applied a zeropoint offset to the COSMOS catalog drawn from a set of realistic offsets Eli provided


resamplings = np.arange(200)

bpzmeans  = {}
truemeans_c = {}

for i in range(5):
    bpzmeans[i] = []
    truemeans_c[i] = []

for r in resamplings:
    d=fitsio.read("%s/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v2.seed"%path + str(r)+"_"+str(r)+"_200.fits")
    dbpz=fitsio.read("%s/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v2.seed"%path + str(r)+"_"+str(r)+"_200.BPZ.fits")
    assert(np.max(np.abs(dbpz['coadd_objects_id']-d['coadd_objects_id']))==0)
    
    # make bins
    bins={}
    bins[0] ,= np.where((dbpz['MEAN_Z']>0.2)&(dbpz['MEAN_Z']<=1.3))
    bins[1] ,= np.where((dbpz['MEAN_Z']>0.2)&(dbpz['MEAN_Z']<=0.43))
    bins[2] ,= np.where((dbpz['MEAN_Z']>0.43)&(dbpz['MEAN_Z']<=0.63))
    bins[3] ,= np.where((dbpz['MEAN_Z']>0.63)&(dbpz['MEAN_Z']<=0.90))
    bins[4] ,= np.where((dbpz['MEAN_Z']>0.90)&(dbpz['MEAN_Z']<=1.3))

    for i in range(5):
        bpzmeans[i].append(np.average(dbpz['Z_MC'][bins[i]]))
        truemeans_c[i].append(np.average(d['redshift'][bins[i]]))



# get some statistics
for i in range(5):
    print("bin",i)
    print("mean BPZ z COSMOS=",np.mean(bpzmeans[i]))
    print("sigma BPZ z COSMOS=",np.std(bpzmeans[i],ddof=1))
    print("mean true z COSMOS=",np.mean(truemeans_c[i]))
    print("SAVE TO TABLE: sigma true z COSMOS=",np.std(truemeans_c[i],ddof=1))
    
for i in range(1,5):
    for j in range(i,5):
        print("bin ",i,j)
        print("(var i)=",np.mean(np.array(truemeans_c[i])**2)-np.mean(truemeans_c[i])**2)
        print("(var j)=",np.mean(np.array(truemeans_c[j])**2)-np.mean(truemeans_c[j])**2)
        print("(covij)=",np.mean(np.array(truemeans_c[i])*np.array(truemeans_c[j]))-np.mean(truemeans_c[j])*np.mean(truemeans_c[i]))

dz_photcalib = np.array([np.std(truemeans_c[i],ddof=1) for i in [1,2,3,4]])
        
print()
print("sigma true z COSMOS goes to https://docs.google.com/document/d/1Bo_zMI1S2F-Han7KkxAS-tDHErKJq0nt9YPu9Vl6EKc")
print("as \"COSMOS photometric calibration uncertainty\"")
print()

def quad(arrays):
  return np.sqrt( sum([np.array(y)**2 for y in arrays]) )


# Boost the priors to reflect new information from clustering
wz_priors= np.array([0.0147, -0.0163, 0.001, 0])
wz_errors = np.array([ 0.026, 0.017, 0.014, 99])
cs_priors = np.array([-0.0119541821982, -0.00952775682453,0.0355976243692,0])
cs_errors = np.array([0.0124704351716,0.0130095999848,0.0112246909963,0])
dz_dis = abs(wz_priors-cs_priors) * (cs_errors/wz_errors)

dz_total = quad([dz_stat[1:], dz_sys, dz_morphology, dz_cv, dz_photcalib, dz_dis])

np.savetxt("%s/cosmos-prior_widths.txt"%args.output,dz_total)

 