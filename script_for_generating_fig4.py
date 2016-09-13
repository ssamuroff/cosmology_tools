import tools.emcee as mc
import tools.plots as pl
reload(pl)
f=pl.fig4()


nf,xmc,ymc=f.extract_line("mcmc", deltaz=True)
ns,xsmc,ysmc=f.extract_line("../fiducial_shear_line/mcmc", deltaz=True)
nr,xrmc,yrmc=f.extract_line("../reduced_line/mcmc", deltaz=True)
tmp,bmc=f.extract_bias_line("mcmc")
tmp,bsmc=f.extract_bias_line("../fiducial_shear_line/mcmc")
tmp,brmc=f.extract_bias_line("../reduced_line/mcmc")

nf,xnrmmc,ynrmmc=f.extract_line("/share/des/disc2/samuroff/mpp/chains/final_postprocessed/robustness/redmagic_deltazprior/halfchain_test/mcmc", deltaz=True)

# From deltaz prior 0.1
nr10,xr10,yr10=f.extract_line("/share/des/disc2/samuroff/mpp/chains/fig4/reduced_line/importance/deltazprior0.10/importance_fromdeltaz0.10/")
xr10,br10=f.extract_bias_line("/share/des/disc2/samuroff/mpp/chains/fig4/reduced_line/importance/deltazprior0.10/importance_fromdeltaz0.10/")
nr6,xr6,yr6=f.extract_line("/share/des/disc2/samuroff/mpp/chains/fig4/reduced_line/importance/deltazprior0.06/importance_fromdeltaz0.06/")
xr6,br6=f.extract_bias_line("/share/des/disc2/samuroff/mpp/chains/fig4/reduced_line/importance/deltazprior0.06/importance_fromdeltaz0.06/")
nr4,xr4,yr4=f.extract_line("/share/des/disc2/samuroff/mpp/chains/fig4/reduced_line/importance/deltazprior0.04/importance_fromdeltaz0.04/")
xr4,br4=f.extract_bias_line("/share/des/disc2/samuroff/mpp/chains/fig4/reduced_line/importance/deltazprior0.04/importance_fromdeltaz0.04/")
nr2,xr2,yr2=f.extract_line("/share/des/disc2/samuroff/mpp/chains/fig4/reduced_line/importance/deltazprior0.02/importance_fromdeltaz0.02/")
xr2,br2=f.extract_bias_line("/share/des/disc2/samuroff/mpp/chains/fig4/reduced_line/importance/deltazprior0.02/importance_fromdeltaz0.02/")
nr1,xr1,yr1=f.extract_line("/share/des/disc2/samuroff/mpp/chains/fig4/reduced_line/importance/deltazprior0.01/importance_fromdeltaz0.01/")
xr1,br1=f.extract_bias_line("/share/des/disc2/samuroff/mpp/chains/fig4/reduced_line/importance/deltazprior0.01/importance_fromdeltaz0.01/")

yr10=np.array([0.0,0.0,0.0,0.0]+list(yr10))
yr6=np.array([0.0]+list(yr6))
br10=np.array([0.0,0.0,0.0,0.0]+list(br10))
br6=np.array([0.0]+list(br6))

yrc = f.composite_line(xr1, [yr1,yr2,yr4,yr6,yr10], [0.01,0.02,0.04,0.06,0.1], yrmc[1:])
brc = f.composite_line(xr1, [br1,br2,br4,br6,br10], [0.01,0.02,0.04,0.06,0.1], brmc[1:])


nr10,xf10,yf10=f.extract_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_line/importance/deltazprior0.10/importance_fromdeltaz0.10/")
xf10,bf10=f.extract_bias_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_line/importance/deltazprior0.10/importance_fromdeltaz0.10/")
nr6,xf6,yf6=f.extract_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_line/importance/deltazprior0.06/importance_fromdeltaz0.06/")
xf6,bf6=f.extract_bias_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_line/importance/deltazprior0.06/importance_fromdeltaz0.06/")
nr4,xf4,yf4=f.extract_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_line/importance/deltazprior0.04/importance_fromdeltaz0.04/")
xf4,bf4=f.extract_bias_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_line/importance/deltazprior0.04/importance_fromdeltaz0.04/")
nr2,xf2,yf2=f.extract_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_line/importance/deltazprior0.02/importance_fromdeltaz0.02/")
xf2,bf2=f.extract_bias_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_line/importance/deltazprior0.02/importance_fromdeltaz0.02/")
nr1,xf1,yf1=f.extract_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_line/importance/deltazprior0.01/importance_fromdeltaz0.01/")
xf1,bf1=f.extract_bias_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_line/importance/deltazprior0.01/importance_fromdeltaz0.01/")

yc = f.composite_line(xf10, [yf1,yf2,yf4,yf6,yf10], [0.01,0.02,0.04,0.06,0.1], ymc[1:])
bc = f.composite_line(xf10, [bf1,bf2,bf4,bf6,bf10], [0.01,0.02,0.04,0.06,0.1], bmc[1:])

nr10,xfs10,yfs10=f.extract_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_shear_line/importance/deltazprior0.10/importance_fromdeltaz0.10/")
xfs10,bfs10=f.extract_bias_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_shear_line/importance/deltazprior0.10/importance_fromdeltaz0.10/")
nr6,xfs6,yfs6=f.extract_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_shear_line/importance/deltazprior0.06/importance_fromdeltaz0.06/")
xfs6,bfs6=f.extract_bias_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_shear_line/importance/deltazprior0.06/importance_fromdeltaz0.06/")
nr4,xfs4,yfs4=f.extract_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_shear_line/importance/deltazprior0.04/importance_fromdeltaz0.04/")
xfs4,bfs4=f.extract_bias_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_shear_line/importance/deltazprior0.04/importance_fromdeltaz0.04/")
nr2,xfs2,yfs2=f.extract_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_shear_line/importance/deltazprior0.02/importance_fromdeltaz0.02/")
xfs2,bfs2=f.extract_bias_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_shear_line/importance/deltazprior0.02/importance_fromdeltaz0.02/")
nr1,xfs1,yfs1=f.extract_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_shear_line/importance/deltazprior0.01/importance_fromdeltaz0.01/")
xfs1,bfs1=f.extract_bias_line("/share/des/disc2/samuroff/mpp/chains/fig4/fiducial_shear_line/importance/deltazprior0.01/importance_fromdeltaz0.01/")

ysc = f.composite_line(xfs10, [yfs1,yfs2,yfs4,yfs6,yfs10], [0.01,0.02,0.04,0.06,0.1], ysmc[1:])
bsc = f.composite_line(xfs10, [bfs1,bfs2,bfs4,bfs6,bfs10], [0.01,0.02,0.04,0.06,0.1], bsmc[1:])





# MCMC points
f.make(xmc,ymc,probe="wl+ggl+lss", extra_label="(fiducial syst.)", linestyle="D",label=True, xlim_lower=0.0, ylim_upper=0.049, lw=3.5)
f.make(xmc,ysmc,probe="wl", extra_label="(fiducial syst.)", linestyle="o",label=True, xlim_lower=0.0, ylim_upper=0.049, lw=3.5, colour="slateblue")
f.make(xmc,yrmc,probe="wl+ggl+lss", extra_label="(reduced syst.)", linestyle="^",label=True, xlim_lower=0.0, ylim_upper=0.049, colour="forestgreen", lw=3.5)
#f.make(xrsmc,yrsmc,probe="wl", extra_label="(reduced syst.)", linestyle="o",label=False, xlim_lower=0.0, ylim_upper=0.044, colour="slateblue", lw=2.5)

f.make(xnrmmc,ynrmmc,probe="wl+ggl+lss", extra_label="(non redmagic)", linestyle=">",label=True, xlim_lower=0.0, ylim_upper=0.049, lw=3.5, colour="darkred")
f.make(xnrmmc,ynrmmc,probe="wl+ggl+lss", extra_label="(non redmagic)", linestyle="-",label=False, xlim_lower=0.0, ylim_upper=0.049, lw=3.5, colour="darkred")

# MCMC points
f.make(xmc,ymc,probe="wl+ggl+lss", extra_label="(fiducial syst.)", linestyle="-",label=False, xlim_lower=0.0, ylim_upper=0.049, lw=3.5)
f.make(xmc,ysmc,probe="wl", extra_label="(fiducial syst.)", linestyle="-",label=False, xlim_lower=0.0, ylim_upper=0.049, lw=3.5, colour="slateblue")
f.make(xmc,yrmc,probe="wl+ggl+lss", extra_label="(reduced syst.)", linestyle="-",label=False, xlim_lower=0.0, ylim_upper=0.049, colour="forestgreen", lw=3.5)
#f.make(xrsmc,yrsmc,probe="wl", extra_label="(reduced syst.)", linestyle="o",label=False, xlim_lower=0.0, ylim_upper=0.044, colour="slateblue", lw=2.5)


#Biases 
f.make(xmc,abs(bsmc),probe="wl", extra_label="(fiducial syst.)", linestyle=":",label=False, xlim_lower=0.0, ylim_upper=0.049, lw=2.5, colour="slateblue")
f.make(xmc,abs(bmc),probe="wl+ggl+lss", extra_label="(fiducial syst.)", linestyle="-.",label=False, xlim_lower=0.0, ylim_upper=0.049, lw=2.5)


#Importance lines
f.make(xf1[2:],yc[2:],probe="wl+ggl+lss", extra_label="(fiducial syst.)", linestyle="-",label=False, xlim_lower=0.0, ylim_upper=0.049, lw=2.5)
f.make(xfs1[2:],ysc[2:],probe="wl", extra_label="(fiducial syst.)", linestyle="-",label=False, xlim_lower=0.0, ylim_upper=0.049, lw=2.5, colour="slateblue")
f.make(xr1[4:],yrc[4:],probe="wl+ggl+lss", extra_label="(reduced syst.)", linestyle="-",label=False, xlim_lower=0.0, ylim_upper=0.049, colour="forestgreen", lw=2.5)

#f.make(xrs10[2:],yrs10[2:],probe="wl", extra_label="(reduced syst.)", linestyle="-",label=False, xlim_lower=0.0, ylim_upper=0.044, colour="slateblue", lw=2.5)

#Some temporary lines pending the 0.02 chains
fid=0.8273733018687
plt.plot([0.0, xf1[2]], np.array([ymc[0], yc[2]])/fid, "-", lw=2.5, color="purple")
plt.plot([0.0, xr1[4]], np.array([yrmc[0], yrc[4]])/fid, "-", lw=2.5, color="forestgreen")
plt.plot([0.0, xfs1[2]], np.array([ysmc[0], ysc[2]])/fid, "-", lw=2.5, color="slateblue")


plt.show()
