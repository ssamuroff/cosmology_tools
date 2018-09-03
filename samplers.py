import numpy as np
import scipy as sp
#import astropy.io.fits as pyfits
import os, pdb
import pylab as plt

labels={"sigma8_input": "$\sigma_8$", "omega_m": "$\Omega _m$", "a": "$A$", "eta": "$\eta$"}

class sampler:
    def __init__(self, fil):
        self.F = np.loadtxt(fil)
        self.param = []
        self.par = {}
        f = open(fil)
        varpar = f.readline().split("\t")
        varied=[]
        likecolumn=False
        for i in varpar:
            if "like" not in i.lower():
                varied.append(i.split("--")[1].replace("\n", ""))
            else:
                likecolumn=True

            # Extract the sampler name
            self.name = f.readline().split("=")[1].replace("\n", "")

        # Read in values config section
        val = f.read().split('VALUES_INI')[1]
        f.close()

        # Separate the lines
        val = val.replace('## ', '').split('\n')

        # Extract the varied parameters
        for c in val:
            if ('=' in c) and (c.split(';')[0].count('.')>1):
                param_name = c.split('=')[0]. replace(' ', '')
                if param_name not in varied:
                    continue
                if param_name in self.par.keys():
                    param_name+='_2'

                self.param += [param_name]
                self.par[param_name] = {}

                param_range = c.split('=')[1].split(';')[0]
                i = 0
                for p in param_range.split(' '):
                    try: 
                        num = float(p)
                        if i==0:
                            self.par[param_name]['min'] = num
                        elif i==1:
                            self.par[param_name]['centre'] = num
                        elif i==2:
                            self.par[param_name]['max'] = num
                        i=i+1
                    except:
                        continue

                print('Found %d x %d results file' %self.F.shape)
                print('Free parameters:', self.param )

                if likecolumn:
                    self.extract_samples(likecolumn)

    def extract_samples(self, likecolumn):
        print("separating columns")
        cols = self.F.T
        for i, p in enumerate(self.param):
            setattr(self, p, cols[i])

            if likecolumn:
                setattr(self, "like", cols[-1])

class emcee(sampler):
    def scatterplt(self, par1=None, par2=None, savedir=None):
        plt.scatter(getattr(self, par1), getattr(self, par2), c=getattr(self, 'like'))
        try:
            plt.xlabel(labels[par1])
        except:
            plt.xlabel(par1)
        try:
            plt.ylabel(labels[par2])
        except:
            plt.ylabel(par2)
        plt.colorbar()
        if not savedir:
            plt.show()
        else:
            plt.savefig(savedir+"sample_points-%s-%s.png"%(par1,par2))


    def find_degeneracy(self, par1, par2, cent=False):
        import scipy.optimize as opt
        x = getattr(self, par1)
        if cent:
            x0 = self.par[par1]["centre"]
        else:
            x0 = 1.
        y = getattr(self, par2)
        par, cov = opt.curve_fit(powerlaw, x/x0, y)
        print(par, cov)

        xfit = np.linspace(x.min(), x.max(), 100)
        yfit = powerlaw(xfit, par[0], par[1])
        plt.plot(xfit,yfit)
        plt.scatter(getattr(self, par1), getattr(self, par2), c=getattr(self, 'like'))
        plt.show()

def powerlaw(x, alpha, c):
    return c*x**(alpha)
