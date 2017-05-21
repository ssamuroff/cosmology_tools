import fitsio as fi
import pylab as plt
import os, yaml, argparse
import glob, copy
import scipy as sp
import numpy as np

description = ''
parser = argparse.ArgumentParser(description=description, add_help=False)
parser.add_argument('--config', type=str, action='store')
args = parser.parse_args()


def load_and_subsample_xis(directory, nbins=np.inf, match_template=False):
    
    if match_template:
        template_fits = fi.FITS("/home/samuroff/hoopoe_paper/cosmology/covariance/2pt_NG.fits")
        theta_target = template_fits["xip"].read()["ANG"]
        bin1 = template_fits["xip"].read()["BIN1"]
        bin2 = template_fits["xip"].read()["BIN2"] 
        template_fits.close()

    theta = np.loadtxt(directory+"/shear_xi/theta.txt") * 180.0 / np.pi * 60.0
    select = (theta>2) & (theta<300)
    if (theta[select].size>nbins):
        thin = int(np.floor(theta[select].size*1./nbins))
        print "subsampling arrays to give %d leave angular bins"%theta.size
    else:
        thin = 1
    xi = {}
    xi["+"] = {}
    xi["theta"] = {}
    xi["-"] = {}
    files = glob.glob(directory+"/shear_xi/xiplus*")
    combinations = [(int(f.split("xiplus_")[1][0]), int(f.split("xiplus_")[1][2]) ) for f in files]
    for (i,j) in combinations:
        print i, j
        xip = np.loadtxt(directory+"/shear_xi/xiplus_%d_%d.txt"%(i,j) )[select][::thin]
        xim = np.loadtxt(directory+"/shear_xi/ximinus_%d_%d.txt"%(i,j) )[select][::thin]
        if match_template:
            t0 = theta_target[(bin1==j) & (bin2==i)]
            interpolate_xip = sp.interpolate.interp1d(np.log(theta[select]), np.log(xip), kind="cubic")
            xip = np.exp(interpolate_xip(np.log(t0)))
            interpolate = sp.interpolate.interp1d(np.log(theta[select]), np.log(xim), kind="cubic")
            xim = np.exp(interpolate(np.log(t0)))
            xi["theta"][(i,j)] = t0


        xi["+"][(i,j)] = xip        
        xi["-"][(i,j)] = xim

    if not match_template:
        x = theta[select][::thin]
    else:
        x = theta_target
    return x, xi


def load_real_xis(info):
    # Load angles, bins, values - should already be in arcminutes
    ptheta, pbin1, pbin2, xip = np.loadtxt(info["xip"]).T
    mtheta, mbin1, mbin2, xim = np.loadtxt(info["xim"]).T

    pbin1+=1
    pbin2+=1
    mbin1+=1
    mbin2+=1

    # should be the same
    if not (ptheta==mtheta).all():
        raise ValueError("Angular binning is inconsistent. Please check the xip/xim files match correctly.")
    else:
        theta = ptheta

    xi = {}
    xi["+"] = {}
    xi["theta"] = {}
    xi["-"] = {}

    for (i,j) in zip(pbin1,pbin2):
        if ((i,j) in xi["+"].keys()):
            continue

        print i, j
        pselect = (pbin1==i) & (pbin2==j) 
        xi["+"][(i,j)] = xip[pselect]
        xi["theta"][(i,j)] = theta[pselect]
        mselect = (mbin1==i) & (mbin2==j)   
        xi["-"][(i,j)] = xim[mselect]

    
    return theta, xi

class xi_template:
    def __init__(self, x, xi):
        print "Setting up xi bias template"
        self.theta = x
        self.correlations = copy.deepcopy(xi)
        print "setting modifier to zero"
        for name in ["+", "-"]:
            for pair in self.correlations[name].keys():
                print name, pair
                self.correlations[name][pair]=np.zeros(self.correlations[name][pair].size)

    def get_mm_correction(self, mmtheta, mmxi, xi, varxi, method="power_law"):
        bf, cov = fit_scale_dependent_m(mmtheta, mmxi, varxi, method)

        print "Generating scale dependent mm correction"
        for name in ["+", "-"]:
            for pair in self.correlations[name].keys():
                print name, pair
                mvec = np.loadtxt("/home/samuroff/shear_pipeline/plot_dump/datavecs/m-vs-z-dvec-fullcat-uncalibrated-desbins.txt").T[-2]
                m2 = np.ones(xi["theta"][pair].size) * mvec[pair[1]-1]
                m1 = np.ones(xi["theta"][pair].size) * mvec[pair[0]-1]

                if (method=="power_law"):
                    dxi_mm = curve(xi["theta"][pair], bf[0], bf[1])
                elif (method=="linear"):
                    dxi_mm = line(xi["theta"][pair][(xi["theta"][pair]<=27)], bf[0], bf[1])
                    dxi_mm = np.hstack((dxi_mm, np.zeros(xi["theta"][pair][(xi["theta"][pair]>27)].size)))
                    
                m1_corr = np.sqrt(m1*m1 + dxi_mm) * m1/abs(m1)
                m2_corr = np.sqrt(m2*m2 + dxi_mm) * m2/abs(m2)
                self.correlations[name][pair] = (1 + m1_corr + m2_corr + m1_corr*m2_corr)/(1 + m1 + m2 + m1*m2)
                import pdb ; pdb.set_trace()
                #dxi_mm * xi[name][pair]

    def apply(self, xi):
        xi_biased = copy.deepcopy(xi)
        print "Applying shift and rescaling"

        for name in ["+", "-"]:
            for pair in self.correlations[name].keys():
                
                xi_biased[name][pair] *= self.correlations[name][pair]
                if self.rescaling!=None:
                    xi_biased[name][pair] *= self.rescaling[name][pair]
                    mr = self.rescaling[name][pair].mean()
                else:
                    mr=1

                print "correlation : ",name, "redshift bin combination : (%d %d)"%pair, "mean rescaling : %3.3f"%mr, "mean shift : %3.3e"%self.correlations[name][pair].mean()
        return xi_biased

def get_neighbours_mean_m_correction(filename, theta, xi):
    data = np.loadtxt(filename).T
    z, m, merr = data[0], data[-2], data[-1]

    print "Applying shift in mean m per bin"

    rescaling = copy.deepcopy(xi)

    for name in ["+", "-"]:
        for (i,j) in xi[name].keys():
            print name, i, j
            m1 = m[i-1]
            m2 = m[j-1]
            rescaling[name][(i,j)] = (1+m1)*(1+m2)

    return rescaling

def fit_scale_dependent_m(theta, xi, varxi, method="power_law"):
    global_mean_m = 0.15126887361063673

    if (method=="power_law"): bf = sp.optimize.curve_fit(curve, theta[2:], xi[2:]-global_mean_m*global_mean_m, p0=[1,0.02])
    if (method=="linear"): bf = sp.optimize.curve_fit(line, theta[2:7], xi[2:7]-global_mean_m*global_mean_m)
    print "best fit curve parameters:" 
    print "a = %f"%bf[0][0]
    print "b = %f"%bf[0][1]

    #import pdb ; pdb.set_trace()

    return bf

def line(x,a,b):
    return np.log10(x)*a+b

def curve(x,a,b):
    return b/(x**a)

def export(filename, theta, xi, nofz, covariance):

    template_fits = fi.FITS("/home/samuroff/hoopoe_paper/cosmology/covariance/2pt_NG.fits")

    print "will export datavector to %s"%filename

    hdus=["xip", "xim", "nz_source"]

    out_xip, out_xim = package_xi(theta, xi)

    out = fi.FITS(filename, "rw")
    out.write(nofz.info)
    out[-1].write_key("EXTNAME", "nz_source")
    for key in nofz.header.keys():
        out[-1].write_key(key, nofz.header[key])

    out.write(out_xip)
    out[-1].write_key("EXTNAME", "xip")
    header = template_fits["xip"].read_header()
    for key in header.keys():
        out[-1].write_key(key, header[key])

    out.write(out_xim)
    out[-1].write_key("EXTNAME", "xim")
    header = template_fits["xim"].read_header()
    for key in header.keys():
        out[-1].write_key(key, header[key])

    out.write(covariance.info)
    out[-1].write_key("EXTNAME", "COVMAT")
    for key in covariance.header.keys():
        out[-1].write_key(key, covariance.header[key])

    
    out.close()

def package_xi(theta, xi):
    dt = [('BIN1', '>i8'), ('BIN2', '>i8'), ('ANGBIN', '>i8'), ('VALUE', '>f8'), ('ANG', '>f8')]
    bin1=[]
    bin2=[]
    angbin=[]
    tbins = np.arange(0,xi["theta"][(1,1)].size, 1) 
    xip=[]
    xim=[]
    ang=[]

    for (i,j) in xi["+"].keys():
        bin1.append([i]*xi["theta"][(i,j)].size)
        bin2.append([j]*xi["theta"][(i,j)].size)
        xip.append(xi["+"][(i,j)])
        xim.append(xi["-"][(i,j)])
        angbin.append(tbins)
        ang.append(xi["theta"][(i,j)])

    out_xip = np.zeros(np.concatenate(bin1).size, dtype=dt)
    out_xim = np.zeros(np.concatenate(bin1).size, dtype=dt)

    out_xip["BIN1"] = np.concatenate(bin1)
    out_xip["BIN2"] = np.concatenate(bin2)
    out_xip["ANG"] = np.concatenate(ang)
    out_xip["ANGBIN"] = np.concatenate(angbin)
    out_xip["VALUE"] = np.concatenate(xip)

    out_xim["BIN1"] = np.concatenate(bin1)
    out_xim["BIN2"] = np.concatenate(bin2)
    out_xim["ANG"] = np.concatenate(ang)
    out_xim["ANGBIN"] = np.concatenate(angbin)
    out_xim["VALUE"] = np.concatenate(xim)
    return out_xip, out_xim





class y1nz:
    def __init__(self, filename, extname="nz_source"):
        print "Loading redshift distributions from %s"%filename
        self.info = fi.FITS(filename)[extname].read()
        self.header = fi.FITS(filename)[extname].read_header()
        self.bins = [name for name in self.info.dtype.names if "BIN" in name]
    def modify(self, info):
        if isinstance(info["truncate"],float):
            z0 = info["truncate"]
            print "WARNING : Will truncate at z=%3.3f"%z0
            select = self.info["Z_MID"]<=z0

            nz = self.info["Z_MID"][select].size

            newinfo = np.zeros(nz, dtype=[('Z_LOW', '>f8'), ('Z_MID', '>f8'), ('Z_HIGH', '>f8'), ('BIN1', '>f8'), ('BIN2', '>f8'), ('BIN3', '>f8'), ('BIN4', '>f8')])

            print "Retaining %d / %d points"%(nz, self.info["Z_MID"].size) 

            for name in self.info.dtype.names:
                newinfo[name] = self.info[name][select]

            self.info = newinfo
            self.header["NAXIS2"] = nz

class y1cov:
    def __init__(self, filename, extname="COVMAT"):
        print "Loading covariance matrix from %s (%s)"%(filename,extname)
        self.info = fi.FITS(filename)[extname].read()[:,:]
        self.header = fi.FITS(filename)[extname].read_header()
        print "Found dimensions %dx%d"%self.info.shape
       
config = yaml.load(open(args.config))

print "Will use %s data"%config["data_type"]
if (config["data_type"]=="mock"):
    theta, xi = load_and_subsample_xis(config["2pt"]["xip"], match_template=True)
else:
    theta, xi = load_real_xis(config["2pt"])


print "Read mm two point function with %d angular bins"%xi["+"][(1,1)].size
print "min scale = %3.3f max scale = %3.3f"%(theta.min(), theta.max())

print "Redshift distributions from %s"%config["nofz"]
nofz = y1nz(config["nofz"])
nofz.modify(config["modifications"]["nofz"])

cov = y1cov(config["covariance"])

export(config["output"], theta, xi, nofz, cov)
