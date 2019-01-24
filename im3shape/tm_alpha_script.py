# Do toy model alpha estimate

from tools.im3shape import basic as i3s
import tools.shapes as s
import numpy as np
m=s.meds_wrapper("/share/des/disc8/cambridge/meds/DES2111+0043-r-sim-ohioA6-meds-y1a1-beta.fits.fz", model="disc")


params = {
    "Rc" :1.67,
    "Rn" :1.17,
    "fc" :1022.70,
    "fn" :555.21,
    "Rp" :1.3,
    "dgn" : 5.}

# m as a fn of PSF ellipticity
m1ofp=[]
m2ofp=[]

angles=np.linspace(0,2*np.pi,61)[:-1]
pshears=np.linspace(-0.08,0.08,9)
gshears=np.array([-0.05,0.05])

for pe1 in pshears:
    e1ofg=[]
    e2ofg=[]
    for g1 in gshears:
        e1vec=[]
        e2vec=[]
        for i,theta in enumerate(angles):  
            x = params["dgn"]*np.cos(theta)
            y = params["dgn"]*np.sin(theta)
            gal,psf=i3s.setup_simple(boxsize=48,shear=(g1,0.0), psf_ellipticity=(pe1,0), psf_size=params["Rp"],  size=params["Rc"], neighbour_ellipticity=(0.0,0.0), neighbour_flux=params["fn"], flux=params["fc"], neighbour_size=params["Rn"], neighbour=[x,y], opt=m.options)
            ra, transform, res, model, img=i3s.i3s([gal],[psf], meds=m, return_all_i3s=True)
            e1vec.append(res.e1)
            e2vec.append(res.e2)
        e1ofg.append([np.mean(e1vec), np.std(e1vec)])
        e2ofg.append([np.mean(e2vec), np.std(e2vec)])
    d1 = np.array(e1ofg).T[0] - gshears
    d2 = np.array(e2ofg).T[0]
    m1ofp.append( (d1[1]-d1[0]) /(shears[1]-shears[0]) )
    m2ofp.append( (d2[1]-d2[0]) /(shears[1]-shears[0]) )


# alpha due to neighbours
# eobs as a fn of psf ellipticit


e1ofp=[]
e2ofp=[]
for pe1 in pshears:
    e1vec=[]
    e2vec=[]
    for i,theta in enumerate(angles):  
        x = params["dgn"]*np.cos(theta)
        y = params["dgn"]*np.sin(theta)
        gal,psf=i3s.setup_simple(boxsize=48,shear=(0.0,0.0), psf_ellipticity=(pe1,0), psf_size=params["Rp"],  size=params["Rc"], neighbour_ellipticity=(0.0,0.0), neighbour_flux=params["fn"], flux=params["fc"], neighbour_size=params["Rn"], neighbour=[x,y], opt=m.options)
        ra, transform, res, model, img=i3s.i3s([gal],[psf], meds=m, return_all_i3s=True)
        e1vec.append(res.e1)
        e2vec.append(res.e2)
    e1ofp.append([np.mean(e1vec), np.std(e1vec)])
    e2ofp.append([np.mean(e2vec), np.std(e2vec)])
    


# Plotting

dm=m1ofp-np.array(m1ofp)[pshears==0][0]

plt.close()
fig, ax = plt.subplots()
ax.plot(pshears,dm*1e3, "o",color="purple", mec="none")
prior_width=0.025
ax.axhline(prior_width*1e3,color="blue", ls="--")
ax.axhline(-prior_width*1e3,color="blue", ls="--")
ax.axhspan(-mbias["m"][1]*1e3,mbias["m"][1]*1e3,color="blue", alpha=0.1)
plt.xlabel("PSF Ellipticity $e^\mathrm{PSF}$")
ax.axhline(0,color="k")
ax.axvline(0,color="k")
plt.xlabel("PSF Ellipticity $e^\mathrm{PSF}$", fontsize=22)
plt.ylabel(r"Excess Multiplicative Bias $\Delta m$ / $\times 10^{-3}$", fontsize=22)
ax.set_xlim(-0.08,0.08)

# left, bottom, width, height
ax2 = fig.add_axes([0.2, 0.16, 0.28, 0.25])
ax2.plot(pshears,dm*1e3, "o",color="purple", mec="none")
ax2.set_ylim(-1.5,0.25)
ax2.set_xlim(-0.08,0.08)
ax2.axhline(0,color="k")
ax2.axvline(0,color="k")
ax2.set_xticks([-0.08,-0.04,0.00,0.04,0.08])
ax2.set_yticks([-1.5, -1, -0.5,0.00])
ax2.tick_params(labelsize=12)


plt.savefig("dm-vs-psfe_neighbours-toy_model.pdf")


   
