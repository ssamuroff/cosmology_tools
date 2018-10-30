import numpy as np
import galsim
import astropy.table as tb
from py3shape import structs
import math
import pylab as plt
plt.switch_backend("agg")
import fitsio
import glob, argparse, os
import tools.shapes as s
import py3shape as p3s
from py3shape import utils
from tools.im3shape import basic as i3s
import tools.diagnostics as di


class mcmc_toy_model:
    def __init__(self):
        print "Toy model to illustrate neighbour bias."

    def get_realisation(self, samples, icent, ineigh, ireal):
        print "-------------------------------------------------"
        print "Generating Realisation %d"%ireal
        print "-------------------------------------------------"

        self.params["fc"] = samples["flux"][icent]
        self.params["fn"] = samples["neighbour_flux"][icent]
        self.params["dgn"] = samples["nearest_neighbour_pixel_dist"][icent]
        self.params["Rc"] = samples["hlr"][icent]
        self.params["Rn"] = samples["neighbour_hlr"][icent]
        self.params["psf_size"] = samples["psf_size"][icent]

        for p in self.params.keys(): print "%s : %3.3f"%(p, self.params[p])

    def sanity_check(self):
        insane = False

        for p in self.params.keys():
            if (self.params[p] < self.priors[p][0]) or (self.params[p] > self.priors[p][1]):
                print "Parameter %s is outside the allowed bounds (%3.3f-%3.3f)"%(p, self.priors[p][0], self.priors[p][1])
                insane = True
        else:
            print 'Parameter values are allowed.'
            insane = False

        return insane

    def run(self, distributions, niterations=2000, filename="mc_toy_model-results", size=1, rank=0):
        # Highest level loop - Sets model parameters
        #-----------------------------------------------------------------------
        
        self.m=[]
        self.centroid=[]

        print "Setting up model"

        meds = s.meds_wrapper("/share/des/disc8/cambridge/meds/DES2111+0043-r-sim-ohioA6-meds-y1a1-beta.fits.fz")


        shears = [-0.02,0.02]
        angles = np.linspace(0,2*np.pi,31)[:-1]

        self.params={"fc":0,"fn":0,"dgn":0,"Rc":0,"Rn":0, "psf_size":0}
        self.priors={"fc":[200,8000],"fn":[2,9000],"dgn":[1,70],"Rc":[0.1,6.0],"Rn":[0.1,6.0], "psf_size":[0.1,4.0]}
        index_c = np.random.choice(distributions.size, niterations*50)
        index_n = np.random.choice(distributions.size, niterations*50)

        idone=0

        for ireal, (icent, ineigh) in enumerate(zip(index_c, index_n)):
            if ireal%size!=rank:
                continue
            self.get_realisation(distributions, icent, ineigh, ireal)
            outside_allowed =  self.sanity_check()
            if outside_allowed:
                continue

            if idone>niterations: continue

            evec_g=[]
            restart=False
            # Second level loop - input shear
            #-----------------------------------------------------------------------
            print "Will evaluate m using %d shear values"%len(shears)
            for ishear, g in enumerate(shears):
                print "g = (%2.2f, 0.00)"%g

                evec_t = []
                centroid=[]

                # Third level loop - neighbour position
                #-----------------------------------------------------------------------
                print "Will use %d neighbour angles"%len(angles)
                for ipos, theta in enumerate(angles):
                    if restart: continue

                    x = self.params["dgn"]*np.cos(theta)
                    y = self.params["dgn"]*np.sin(theta)
                    print "theta = %2.3f degrees, position = (%3.2f,%3.2f)"%(theta*60., x, y)

                    gal,psf = i3s.setup_simple(boxsize=96,shear=(g,0.0), psf_size=self.params["psf_size"],  size=self.params["Rc"], neighbour_ellipticity=(0.0,0.0), neighbour_flux=self.params["fn"], flux=self.params["fc"], neighbour_size=self.params["Rn"], neighbour=[x,y], opt=meds.options)
                    res = i3s.i3s([gal.array],[psf], meds=meds)

                    evec_t.append([res.e1, res.e2])
                    centroid.append(np.sqrt(res.ra_as*res.ra_as + res.dec_as*res.dec_as))
                

                meane1 = np.array(evec_t).T[0].mean()
                meane2 = np.array(evec_t).T[1].mean()
                evec_g.append([meane1, meane2])
            

            # Finally we have a vector, containing one mean measured shape for each of the input shear values
            # Calculate m
            residual_e1 = np.array(evec_g).T[0] - np.array(shears)
            residual_e2 = np.array(evec_g).T[1]

            m=(residual_e1[-1]-residual_e1[0])/(shears[-1]-shears[0])


            print "---------------------- m=%f"%m
            print centroid

            self.m.append([(residual_e1[-1]-residual_e1[0])/(shears[-1]-shears[0]), (residual_e2[-1]-residual_e2[0])/(shears[-1]-shears[0])]) 
            self.centroid.append([np.array(centroid).mean(), np.array(centroid).max()])
            if abs(self.m[-1][0])>2: continue
            if (self.m[-1][0]<-0.01) and (np.array(centroid).mean()<1): import pdb ; pdb.set_trace()
            self.write_output_line(filename)
            idone+=1
        
        print "Done all loops"


    def write_output_line(self, filename):
        f = open(filename, "a")
        line = "%3.4f %3.4f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %3.6f %3.6f \n"%(self.centroid[-1][0], self.centroid[-1][1], self.params["Rc"], self.params["Rn"], self.params["fc"], self.params["fn"], self.params["dgn"], self.params["psf_size"], self.m[-1][0], self.m[-1][1])
        f.write(line)
        f.close()

    def do_position_loop(self, proxy):
        """Loop over all possible angular positions for this neighbour realisation."""
        ang=[]
        m_theta={}
        theta= np.linspace(0,2,10)*np.pi
        snr = []

        for i,t in enumerate(theta):
            print "  --  Level 2 iteration: %d (theta=%f)"%(i,t)

            m_theta[1]=[]
            m_theta[2]=[]

            x = self.dneigh*np.cos(t)
            y = self.dneigh*np.sin(t)

            sn, m1,m2 = self.simple_m(proxy, x, y)
            print "At SNR = %4.1f got m = (%2.2f,%2.2f)"%(sn,m1,m2)
            m_theta[1].append(m1)
            m_theta[2].append(m2)
            ang.append(t)
            snr.append(sn)

        # Compute the mean over all positions
        self.m[1].append(np.mean(m_theta[1]))
        self.m[2].append(np.mean(m_theta[2]))
        self.snr.append(np.mean(snr))

    def simple_m(self, proxy, x, y):

        y_g = []
        snr = []
        shear = [-0.02,0.02]
        for g in shear:
            gal, psf = setup_simple(boxsize=48, shear=(self.e1+g,self.e2), flux=proxy, psf_ellipticity=(0.0,0.0), psf_size=self.psf_size, neighbour_ellipticity=(0.0,0.0), neighbour=[x,y], neighbour_flux=self.neighbour_flux, neighbour_size=self.neighbour_size, wcs=self.wcs, opt=self.opt)
            res, cat = run(gal, psf, opt=self.opt, return_cat=True, show=False)
            snr.append(cat.snr)
            y_g.append(res.e1 - self.e1 - g)

        m1 = (y_g[1]-y_g[0])/(shear[1]-shear[0])
        snr = np.mean(snr)

        y_g = []
        for g in shear: 
            gal, psf = setup_simple(boxsize=48, shear=(self.e1,self.e2+g), flux=proxy, psf_ellipticity=(0.0,0.0), psf_size=self.psf_size, neighbour_ellipticity=(0.0,0.0), neighbour=[x,y], neighbour_flux=self.neighbour_flux, neighbour_size=self.neighbour_size, wcs=self.wcs, opt=self.opt)
            res, cat = run(gal, psf, opt=self.opt, return_cat=True)
            y_g.append(res.e2 - self.e2 -  g)

        m2 = (y_g[1]-y_g[0])/(shear[1]-shear[0])

        return snr, m1, m2

    def show(self, xmin=-0.2,xmax=0.2,noneigh=True):

        fig = plt.figure()
        plt.plot(self.binning,self.g["e1"]*1e4,color="purple",lw=2.5, ls="-",label="with neighbour")
        if noneigh:
            plt.plot(self.binning,self.g["e2"]*1e4,color="steelblue",lw=2.5, ls=":",label="without neighbour")
        plt.axhline(0.0, color="k", linestyle="--", lw=2.5)
        plt.axvline(0.0, color="k", linestyle="--", lw=2.5)
        plt.xlim(xmin,xmax)

        plt.xlabel(r"PSF ellipticity $ e_1 ^{psf}$")
        plt.ylabel(r"mean ellipticity $<e_1>_{\theta} \times 10^4$")
        if noneigh:
            plt.legend(loc="upper left")
        plt.show()

    def alpha(self, component=1):
        sel = abs(self.binning)<0.1
        linear_e = self.g["e%d"%component][sel]
        linear_epsf = self.binning[sel]
        return (linear_e[-1]-linear_e[0])/(linear_epsf[-1]-linear_epsf[0])

    def sample_alpha(self, param_name, samples):
        alpha = []
        print "Will sample alpha at points in %s"%param_name
        for q in samples:
            print q,
            exec "self.generate(%s=%f)"%(param_name,q)
            alpha.append(self.alpha())
            print alpha[-1]

        return samples, alpha



class results:
    def __init__(self, directory="/home/samuroff/hoopoe_paper/toy_model_data", impose_priors=True):
        self.samples=[]
        files = glob.glob("%s/mc_toy_model-results*.txt"%directory)
        for f in files:
            print os.path.basename(f)
            self.samples.append(tb.Table.read(f, format="ascii", names=['Rc','Rn', 'fc', 'fn','dgn','psf_size', 'm', 'm12']))

        self.samples = np.concatenate(self.samples)

        self.priors={"fc":[200,8000],"fn":[0.,12000],"dgn":[1,12],"Rc":[0.1,8.0],"Rn":[0.1,8.0], "psf_size":[0.1,4.0]}
        if impose_priors:
            self.impose_priors()

    def reset_priors(self, priors, name=None):
        if isinstance(priors, dict):
            self.priors = priors
        elif isinstance(priors, dict) and (name is not None):
            self.priors[name] = priors
        else:
            print "Priors given in unrecognised format"

    def impose_priors(self):
        nsamp = self.samples.size
        for q in self.priors.keys():
            select = (self.samples[q]>self.priors[q][0]) & (self.samples[q]<self.priors[q][1])
            ncut = nsamp - len(self.samples[select])
            print "Prior on %s [%3.2f-%3.2f] removes %d (%3.3f percent) samples"%(q, self.priors[q][0], self.priors[q][1], ncut, ncut*100./nsamp )
            self.samples = self.samples[select]

    def generate_binned_data(self, nbin=6):
        self.vec={}
        for name in self.samples.dtype.names: 
            if name in ["m", "m12"]: continue
            print name
            bins=di.find_bin_edges(self.samples[name], nbin)
            tmp=[]
            x=[]           
            for b in zip(bins[:-1], bins[1:]):
                print b
                selection = (self.samples[name]<b[1]) & (self.samples[name]>b[0]) 
                tmp.append([self.samples["m"][selection].mean(), self.samples["m"][selection].std(), self.samples["m"][selection].size])
                x.append(self.samples[name][selection].mean())
            self.vec[name]=(x,tmp)


