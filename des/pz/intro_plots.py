import pylab as plt
import fitsio
import numpy as np
plt.switch_backend("pdf")
plt.style.use("y1a1")

# (1) make some plots purely for orientation
# (1a) plot p(z) of 200 COSMOSes resampled/reweighted with different photometric calibration offsets

for i in range(200):
    # 200 resamplings/reweightings of COSMOS metacal galaxies with 200 different random photometric offsets applied to the COSMOS data
    # each of these catalogs is based on 200.000 galaxies in the metacalibration sample
    # to each one of them we find a (color+magnitude) matching galaxy in the COSMOS sample (some of the COSMOS galaxies appear multiple times)
    # from the galaxy in the actual metacalibration sample we know its shape measurement success / weight
    # from the COSMOS galaxy we know the 'true' redshift
    d=fitsio.read("/share/des/disc4/samuroff/ia/pz_calibration/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v2.seed"+str(i)+"_"+str(i)+"_200.fits")
    plt.hist(d['redshift'],weights=(d['R11']+d['R22']),bins=400,range=(0,2),histtype='step',alpha=0.01,color='black',normed=False)
    # d['redshift'] is the 'true' (30-band photo-z) redshift of galaxies in COSMOS
    # the effective weight of a galaxy in the shape catalog is proportional to R11+R22, the shear response
    
plt.xlabel(r"$z_{\rm true}$")
plt.ylabel(r"$p(z)$")
plt.title("COSMOS matched p(z) with photo-calibration noise")
plt.savefig("cosmos_matched_pz.pdf")
plt.close()


# (1b) checking how frequently galaxies get reused in COSMOS resampling
for i in range(200):
    d=fitsio.read("/share/des/disc4/samuroff/ia/pz_calibration/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v2.seed"+str(i)+"_"+str(i)+"_200.fits")
    unique, idx, counts = np.unique(d['matched_coadd_objects_id'], return_index=True, return_counts=True)
    # matched_coadd_objects_id is a unique ID of the COSMOS galaxy that was used as a match
    # it appears multiple times if a COSMOS galaxy is the best match for multiple metacalibration sample galaxies
    plt.hist(counts,bins=200,range=(0,200),histtype='step',alpha=0.02,color='purple',log=True)

plt.xlabel("repetitions of COSMOS galaxy")
plt.title("a few galaxies are too popular, but it's not too bad")
plt.savefig("cosmos-repetition.pdf")
plt.close()


# (1c) look at the distribution of <z>, in samples with different photometric calibration
m=[]
for i in range(200):
    d=fitsio.read("/share/des/disc4/samuroff/ia/pz_calibration/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v2.seed"+str(i)+"_"+str(i)+"_200.fits")
    m.append(np.average(d['redshift'],weights=(d['R11']+d['R22']))) 
    # weighted mean of true redshift, weight is metacal response of associated metacal sample galaxy
plt.hist(m,bins=20,alpha=0.2)
print("mean z=",np.average(m),"+-",np.std(m))
plt.xlabel(r"$z_{\rm mean}$")
plt.title("metacal-weighted mean z of full sample, 200 photo-calib realiations")
plt.savefig("metacal-weighted-mean-z-fullsamp-200photocalib-realiations.pdf")
plt.close()

# (1d) in one resampling, look at the distribution of mof magnitudes
# note that these can be below detection threshold if a lot of flux gets assigned to a blended neighbor
plt.hist(d['mof_mag_i'],bins=100,log='True',alpha=0.2)
plt.xlim((16,28.))
plt.xlabel("MOF mag i")
plt.savefig("mof-mag-i.pdf")
# I've checked that the highly repeated COSMOS galaxies are bright galaxies that are just rare in COSMOS
plt.close()


# (1e) do the galaxies with high multiplicity matter for the overall p(z)? 
# plot their z histogram only
for i in range(200):
    d=fitsio.read("/share/des/disc4/samuroff/ia/pz_calibration/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v2.seed"+str(i)+"_"+str(i)+"_200.fits")
    unique, idx, counts = np.unique(d['matched_coadd_objects_id'], return_index=True, return_counts=True)
    idx = idx[counts>20]
    plt.hist(d['redshift'][idx],weights=(d['R11'][idx]+d['R22'][idx]),bins=400,range=(0,2),histtype='step',alpha=0.01,color='black',normed=False)

plt.xlabel(r"$z_{\rm true}$")
plt.ylabel(r"$p(z)$")
plt.title("COSMOS matched $p(z)$, multiplicity$>20$")
plt.ylim(0,20)
plt.savefig("repeated-cosmos-galaxies-hist.pdf")
plt.close()
# note the y axis  scale compared to figure (1a) -- this is a non-issue


# (1f) dig a little deeper into the spikes in 1a and see whether we find something odd
# let's look at one of the spikes
d=fitsio.read("output_collated/wl_class.METACAL_MOF.rescaled.slr.cosmos.v2.seed0_0_200.fits")
plt.hist(d['redshift'][(d['redshift']>0.23)&(d['redshift']<0.24)],bins=400,range=(0,2),histtype='step',color='black',normed=False)
plt.savefig("spike-examp.pdf")
plt.close()

e=d[(d['redshift']>0.23)&(d['redshift']<0.24)] # there is a spike between 0.23 and 0.24, select just those objects

plt.hist(e['mof_mag_i'],bins=100,range=(17,25))
plt.xlabel("i")
plt.savefig("modf-mag-hist-spike.pdf")

# these objects have a more or less sane magnitude distribution



# plot their positions in the actual COSMOS field
f=fitsio.read("COSMOSY1V103.Y1A1_METACAL_MOF_D04.v2.fits")
ra=[]
dec=[]
for i in e:
    fg=f[f["coadd_objects_id"]==i["matched_coadd_objects_id"]]
    ra.append(fg["ra"][0])
    dec.append(fg["dec"][0])

plt.scatter(ra,dec,marker="x",alpha=0.3,color="steelblue")
plt.xlabel("ra")
plt.ylabel("dec")
plt.savefig("cosmos-positions.pdf")
# they look kind of clustered, but they're not all in one place, so fine

