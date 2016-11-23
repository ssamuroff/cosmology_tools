from tools import fisher as f

shear=f.fisher("shear-only-fisher-syst+cos.txt")
sg=f.fisher("shear+ggl-fisher-1dshift-syst+cos.txt")
sp=f.fisher("shear+pos-fisher-1dshift-cut-sys+cos.txt")
sgp=f.fisher("shear+ggl+pos-fisher-1dshift-syst+cos.txt")

shear.reset_priors("delta-priors-shear.ini")
sg.reset_priors("delta-priors-ggl.ini")
sgp.reset_priors("delta-priors-ggl.ini")
sp.reset_priors("delta-priors-ggl.ini")

x2,y = f.photoz_prior_width_plot(shear, lens=False)
x2,ysg = f.photoz_prior_width_plot(sg, lens=False)
x2,ysp = f.photoz_prior_width_plot(sp, lens=False)
x2,ysgp = f.photoz_prior_width_plot(sgp, lens=False, show=True)

x,y = f.photoz_prior_width_plot(shear, lens=True)
x,ysg = f.photoz_prior_width_plot(sg, lens=True)
x,ysp = f.photoz_prior_width_plot(sp, lens=True)
x,ysgp = f.photoz_prior_width_plot(sgp, lens=True)

plt.plot(x,y/0.82, 'm-', label='WL')
plt.plot(x,ysg/0.82, 'g-', label='WL+ggl')
plt.plot(x,ysp/0.82, 'r-', label='WL+LSS')
plt.plot(x,ysgp/0.82, 'b-', label='WL+ggl+LSS')
plt.xlim(0.,0.06)
plt.legend(loc='upper left')
plt.xlabel("prior width $\Delta \delta z$")
plt.ylabel("fractional error $\Delta \sigma_8 / \sigma_8$")
plt.show()

plt.plot(x,y/y[0], 'm-', label='WL')
plt.plot(x,ysg/ysg[0], 'g-', label='WL+ggl')
plt.plot(x,ysp/ysp[0], 'r-', label='WL+LSS')
plt.plot(x,ysgp/ysgp[0], 'b-', label='WL+ggl+LSS')
plt.xlim(0.,0.06)
plt.legend(loc='upper left')
plt.xlabel("prior width $\Delta \delta z$")
plt.ylabel("fractional degradation $\Delta \sigma_8 / \Delta \sigma_{8,0}$")
plt.show()

