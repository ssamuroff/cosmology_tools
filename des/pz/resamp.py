import numpy as np
import fitsio as fi
import os,yaml,glob
import argparse
from string import Template as Tmp
import socket
import pylab as plt
plt.switch_backend("pdf")
plt.style.use("y1a1")

def wavg(x,wt):
    return sum(x*wt)/sum(wt)

class resampled_cosmos:
    def __init__(self):
        self.hostname = socket.gethostname()
        if 'fornax' in self.hostname:
            self.path="/share/des/disc4/samuroff/ia/pz_calibration/"
        else:
            self.path="/global/cscratch1/sd/sws/y1/photoz/"

    def load(self):
        f="%s/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.fits"%self.path
        self.rphot=fi.read(f) # resampled photometry, 200000 random objects from the science sample                                   
        self.rmofbpz=fi.read("%s/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.BPZMOF.fits"%self.path) # BPZ run on that MOF photometry  
        self.rmcalbpz=fi.read("%s/cats/wl_class.METACAL_MOF.rescaled.slr.cosmos.v4.BPZ.fits"%self.path) # BPZ run on that MCAL photometry
        self.truebpz = fi.read("%s/hst_fluxes-COSMOSY1V103.Y1A1_COADD_D04.v2.BPZ.fits"%self.path)
    def get_bootstrap_errors(self, niter, t_des, t_cos, R, rmcalbpz,t_bounds):
        red_realisations = []
        blue_realisations = []
        for i in xrange(niter):
            subset = np.random.choice(R.size, R.size/2, replace=False)
            des_not_cosmos = (t_des[subset]<t_bounds[0]) & (t_cos[subset]>t_bounds[1])
            cosmos_not_des = (t_cos[subset]<t_bounds[0]) & (t_des[subset]>t_bounds[1])
            R_red = sum(R[subset][t_des[subset]<t_bounds[0]])
            blue_leakage = sum(R[subset][des_not_cosmos])/R_red
            red_leakage  = sum(R[subset][cosmos_not_des])/R_red
            red_realisations.append(red_leakage)
            blue_realisations.append(blue_leakage)
        return np.std(red_realisations), np.std(blue_realisations)
    def assess_purity(self, bootstrap_iterations=0, t_bounds=[1,1]):
        """We're trying to assess how strongly a sample defined by a cut on observed T_BPZ is affected 
           by leakage
           (ie how much it matters that there are galaxies that would be classified as blue if
           we had HST-quality photometry, but are mis-classified as red because we're using
           noisy, potentially biased observations).
           Use a BPZ run on the HST resampled catalogue to try to gauge this. """
        bounds = [(0.20,1.30),(0.20,0.43),(0.43,0.63),(0.63,0.90),(0.90,1.30)]
        lmask = [ np.where((self.rmcalbpz['MEAN_Z']>a)&(self.rmcalbpz['MEAN_Z']<b))[0] for (a,b) in bounds]
        num = []
        if bootstrap_iterations>0:
            print "Will compute bootstrap errorbars (%d iterations)"%bootstrap_iterations
            print "this may take a little while --  set 'bootstrap_iterations=0' to skip."
        er, eb = 0,0 
        print "Colour bins : T<%3.3f, T>%3.3f"%tuple(t_bounds)
        for j,l in enumerate(lmask):
            print j
            indices = np.array([np.argwhere(self.truebpz["coadd_objects_id".upper()]==i)[0,0] for i in self.rmcalbpz['matched_coadd_objects_id'][l]])
            if not (self.truebpz["coadd_objects_id".upper()][indices]==self.rmcalbpz['matched_coadd_objects_id'][l]).all(): import pdb ; pdb.set_trace()
            t_des = self.rmcalbpz['TEMPLATE_TYPE'][l]
            t_cos = self.truebpz["template_type".upper()][indices]
            R = (self.rphot["R11"]+self.rphot["R22"])[l]
            zr = wavg(self.rmcalbpz["redshift"][l],R)
            R_red = sum(R[t_des<1])
            des_not_cosmos = (t_des<t_bounds[0]) & (t_cos>t_bounds[1])
            cosmos_not_des = (t_cos<t_bounds[0]) & (t_des>t_bounds[1])
            blue_leakage = sum(R[des_not_cosmos])/R_red  # Weight of galaxies which are actually blue according to the BPZ run on the HST fluxes,
                                                                       # but the run on DES thinks are red
            red_leakage  = sum(R[cosmos_not_des])/R_red  # Weight of galaxies which are actually red according to the BPZ run on the HST fluxes,
                                                                       # but the run on DES thinks are blue
            if bootstrap_iterations>0: er, eb = self.get_bootstrap_errors(bootstrap_iterations,t_des,t_cos,R,self.rmcalbpz[l],t_bounds) 
            num.append(np.array([z,blue_leakage, eb, red_leakage, er, R[des_not_cosmos].size, R[cosmos_not_des].size ]) )
        num = np.array(num)
        print num
        return num
    def select(self,gtype):
        cuts = yaml.load(open("%s/cuts/%s.yaml"%(self.path,gtype)))
        self.mcal_mask = np.ones(self.rmcalbpz.size).astype(bool)
        self.mof_mask = np.ones(self.rmofbpz.size).astype(bool)
        for col in cuts.keys():
            print "Cutting %s"%col
            lower,upper,dt = cuts[col].split()
            exec("lower=%s(lower)"%dt)
            exec("upper=%s(upper)"%dt)
            self.mcal_mask = self.mcal_mask & (self.rmcalbpz[col]>lower) & (self.rmcalbpz[col]<upper)
            self.mof_mask = self.mof_mask & (self.rmofbpz[col]>lower) & (self.rmofbpz[col]<upper)


def purity_vs_tlow():
    cosmos = resampled_cosmos()
    cosmos.load()
    t0 = np.linspace(0.3,1,8)
    vec=[]
    for t_low in t0:
        fe = cosmos.assess_purity(bootstrap_iterations=50, t_bounds=[t_low,1])
        vec.append([fe.T[1],fe.T[2],fe.T[3],fe.T[4]])

    return t0, np.array(vec)
       
