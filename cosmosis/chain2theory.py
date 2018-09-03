import numpy as np
import scipy, glob, argparse, yaml
import fitsio as fi
import pylab as plt
import matplotlib
import os
from tools.cosmosis import cosmosis_tools as ct
import tools.emcee as mc
from string import Template as Tm

positions={
    (1,1,"+"):4,
    (1,1,"-"):12,
    (2,1,"+"):3,
    (2,1,"-"):18,
    (3,1,"+"):2,
    (3,1,"-"):24,
    (4,1,"+"):1,
    (4,1,"-"):30,
    (2,2,"+"):9,
    (2,2,"-"):17,
    (3,2,"+"):8,
    (3,2,"-"):23,
    (4,2,"+"):7,
    (4,2,"-"):29,
    (3,3,"+"):14,
    (3,3,"-"):22,
    (4,3,"+"):13,
    (4,3,"-"):28,
    (4,4,"+"):19,
    (4,4,"-"):27}


def get_theory_spectra(i,j,filename):
    theta = np.loadtxt("%s/shear_xi/theta.txt"%filename)
    theta = theta * 60. * 180. / np.pi
    xip = np.loadtxt("%s/shear_xi/xiplus_%d_%d.txt"%(filename,i+1,j+1))
    xim = np.loadtxt("%s/shear_xi/ximinus_%d_%d.txt"%(filename,i+1,j+1))
    return theta, xip, xim

def get_real_spectra(i,j,fits,error=True):
    if fits is None:
        return [],[],[]
    else:
        selectp = (fits["xip"]["BIN1"][:]==j+1) & (fits["xip"]["BIN2"][:]==i+1)
        selectm = (fits["xim"]["BIN1"][:]==j+1) & (fits["xim"]["BIN2"][:]==i+1)

        xp = fits["xip"]["ANG"][:][selectp]
        xm = fits["xim"]["ANG"][:][selectm]

        xip = fits["xip"]["VALUE"][:][selectp]
        xim = fits["xim"]["VALUE"][:][selectm]

        if error:
            cov = fits["COVMAT"][:,:]
            startp,endp = fits["covmat"].read_header()["STRT_0"], fits["covmat"].read_header()["STRT_1"]
            dx = endp - startp
            startm = endp
            endm = startm + dx

            nx = xp.size

            covp = cov[startp:endp,startp:endp]
            covm = cov[startm:endm,startm:endm]

            i0 = nx*(i + j)
            errp = np.diag(covp[i0:(i0+nx),i0:(i0+nx)])
            errp = np.sqrt(errp)
            errm = np.diag(covm[i0:(i0+nx),i0:(i0+nx)])
            errm = np.sqrt(errm)
            if len(errp)==0:
                import pdb ; pdb.set_trace()

        else:
            errp = None
            errm = None

        return (xp,xip,errp), (xm,xim,errm)

def get_modifier(x,dxi,target):
    if len(dxi)==0: return 0
    import scipy.interpolate as interp
    interpolator = interp.interp1d(np.log10(x),dxi) 
    return interpolator(np.log10(target))
    

def cornerplot(theory, data1, data2, dz=[]):

    if len(dz)>0:
        print "will correct photo-z error for corner plot"
        correct=True
    else:
        correct=False
        dz=[0,0,0,0]

    plt.switch_backend("pdf")
    plt.style.use("y1a1")
    matplotlib.rcParams["ytick.minor.visible"]=False
    matplotlib.rcParams["ytick.minor.width"]=0.1
    ni,nj=np.genfromtxt("%s/shear_xi/values.txt"%theory).T[2]
    ni = int(ni)
    nj = int(nj)
    if data1 is not None:
        data1 = fi.FITS(data1)
    if data2 is not None:
        data2 = fi.FITS(data2)

    rows, cols = ni+1, nj+2

    count = 0

    for i in xrange(ni):
        for j in xrange(nj):
            count+=1
            if j>i:
                continue

            print i,j
            xt,xip,xim = get_theory_spectra(i,j,theory)
            xip_a,xim_a = get_real_spectra(i,j,data1)
            xip_b,xim_b = get_real_spectra(i,j,data2)

            posp = positions[(i+1,j+1,"+")]
            ax = plt.subplot(rows,cols,posp)
            ax.annotate("(%d, %d)"%(i+1,j+1), (30,8e-5), textcoords='data', fontsize=12)
            ax.yaxis.set_tick_params(which='minor', left='off', right='off')
            plt.ylim(1e-9,1e-3)
            plt.yscale("log")
            plt.xscale("log")
            if posp==19:
                plt.ylabel(r"$\xi_+(\theta)$", fontsize=12)
                plt.xlabel(r"$\theta$ / arcmin", fontsize=12)
                plt.yticks(fontsize=12)
            else:
                plt.yticks(visible=False)
            plt.xlim(2.2,270)
            plt.xticks([10,100],["10", "100"])
            plt.plot(xt,xip,color="k")

            positive = (xip_a[1]>0)
            #import pdb ; pdb.set_trace()
            dy = get_modifier(dz[0],dz[1][i+1,j+1],xip_a[0][positive])
            plt.errorbar(xip_a[0][positive], xip_a[1][positive]+dy, yerr=xip_a[2][positive], ls="none", marker=".", ecolor="red", markeredgecolor="red", markerfacecolor="none")
            positive = (xip_b[1]>0)
            dy = get_modifier(dz[0],dz[2][i+1,j+1],xip_b[0][positive])
            plt.errorbar(xip_b[0][positive], xip_b[1][positive]+dy, yerr=xip_b[2][positive], ls="none", marker="x", ecolor="steelblue", markeredgecolor="steelblue", markerfacecolor="steelblue")

            posm = positions[(i+1,j+1,"-")]
            ax = plt.subplot(rows,cols,posm)
            ax.annotate("(%d, %d)"%(i+1,j+1), (30,8e-5), textcoords='data', fontsize=12)
            ax.yaxis.set_tick_params(which='minor', left='off', right='off')
            plt.ylim(1e-9,1e-3)
            if posm==30:
                ax.yaxis.set_label_position("right")
                ax.yaxis.set_ticks_position("right")
                plt.yticks(fontsize=12)
                plt.ylabel(r"$\xi_-(\theta)$", fontsize=12)
            else:
                plt.yticks(visible=False)
            plt.yscale("log")
            plt.xscale("log")
            plt.xlim(2.2,270)
            plt.xticks([10,100],["10", "100"])
            plt.plot(xt,xim,color="k")
            positive = (xim_a[1]>0)
            dy = get_modifier(dz[0],dz[3][i+1,j+1],xip_a[0][positive])

            plt.errorbar(xim_a[0][positive], xim_a[1][positive]+dy, yerr=xim_a[2][positive], ls="none", marker=".", ecolor="red", markeredgecolor="red", markerfacecolor="none")
            positive = (xim_b[1]>0)
            dy = get_modifier(dz[0],dz[4][i+1,j+1],xip_b[0][positive])
            plt.errorbar(xim_b[0][positive], xim_b[1][positive]+dy, yerr=xim_b[2][positive], ls="none", marker="x", ecolor="steelblue", markeredgecolor="steelblue", markerfacecolor="steelblue")

    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig("/home/samuroff/tmp3.pdf")


def preprocess(filename, frac=1.0):
    os.system("mkdir -p output/")
    c = mc.chain(filename)
    c.add_column("s8", values="sigma_8*((omega_m/0.3)**2)")
    nsamp = len(c.samples)   
    nburn = int( (1.0 - frac) * nsamp)
    if (nburn!=0):
        c.burn(nburn)
    c.write_columns("output/chain.txt")

def parse_values(filename, args, blind=False):
    mean_vals = np.genfromtxt("%s/means.txt"%filename, dtype=[("name", "S100"),("mean", float), ("std", float)])

    params = {
    "cosmological_parameters":
        {
        "omega_m" :  0.295,
        "h0"       :  0.6881,
        "omega_b"  :  0.0468,
        "n_s"        :  0.9676,
        "A_s"        :  2.260574e-9,
        "omnuh2"     :  0.0006,
        "w"         : -1.0,
        "massive_nu" :  1,
        "massless_nu":  2.046,
        "omega_k"    :  0.0,
        "tau"        :  0.08,
        "wa"         :  0.0},
    "shear_calibration_parameters":
        {
        "m1" : 0.0,
        "m2" : 0.0,
        "m3" : 0.0,
        "m4" : 0.0},
    "intrinsic_alignment_parameters":
        {
        "a"     : 0.0,
        "alpha" : 0.0,
        "z0"    : 0.62},
    "wl_photoz_errors":
        {
        "bias_1" : 0.0,
        "bias_2" : 0.0,
        "bias_3" : 0.0,
        "bias_4" : 0.0}}

    for i, name in enumerate(mean_vals["name"]):
        if (name=="post"):
            continue
        section, param = name.split("--")
        print section, param,
        value = mean_vals["mean"][i]
        params[section][param] = value
        if blind:
            print "XXX"
        else:
            print value

    if args.noia:
        params["intrinsic_alignment_parameters"]["a"]=0.0

    return params

def export_values(mean_values_dict):
    vals = Tm(open("/home/samuroff/local/python/lib/python2.7/site-packages/tools/cosmosis/values_template").read())
    vals_all = {}
    for name1 in mean_values_dict.keys():
        for  name2 in mean_values_dict[name1].keys():
            vals_all[name2] = mean_values_dict[name1][name2]

    values_txt = vals.substitute(vals_all)

    outfile = open("output/values.ini", "wa")
    outfile.write(values_txt)
    outfile.close()

def replace_dz_values(source,values):
    new={}
    for key in values.keys():
        new[key] = {}
        for key2 in values[key].keys():
            new[key][key2] = values[key][key2]
    if source=="none":
        for i in [1,2,3,4]: 
            new["wl_photoz_errors"]["bias_%d"%i]=0
    else:
        priors = "/home/samuroff/priors/%s-v5/mcal_priors.ini"%source
        mean_vals = {'early': [ -0.0218003622646,  -0.0398050264496, -0.008459663380,  -0.0442398688684],
                      'late': [-0.00348994944362, -0.00680842770299, 0.0297324014146, -0.00968806318733]}
        for i in [1,2,3,4]: 
            print mean_vals[source][i-1]
            new["wl_photoz_errors"]["bias_%d"%i] = mean_vals[source][i-1]
    return new

def get_shift():
    dzp = {}
    dzm = {}
    for i in [1,2,3,4]:
        for j in [1,2,3,4]:
            if i<j:
                continue
            x = np.loadtxt("test_run1/shear_xi/theta.txt") * 180. / np.pi * 60.
            xip = np.loadtxt("test_run1/shear_xi/xiplus_%d_%d.txt"%(i,j)).T 
            xip_ref = np.loadtxt("test_run2/shear_xi/xiplus_%d_%d.txt"%(i,j)).T 
            xim = np.loadtxt("test_run1/shear_xi/ximinus_%d_%d.txt"%(i,j)).T 
            xim_ref = np.loadtxt("test_run2/shear_xi/ximinus_%d_%d.txt"%(i,j)).T
            dzp[(i,j)] = xip - xip_ref
            dzm[(i,j)] = xim - xim_ref

    return x,dzp,dzm


def main(args):
    # postprocess the specified chain
    preprocess(args.chain, frac=args.fkeep)
    ct.postprocess("output/chain.txt", "output")

    # Read and convert the mean values to an ini file cosmosis can run with
    mean_values_dict = parse_values("output", args, blind=args.blind)
    export_values(mean_values_dict)

    # Finally call cosmosis
    ini = "/home/samuroff/local/python/lib/python2.7/site-packages/tools/cosmosis/ini/des-test.ini"
    ct.run_cosmosis(ini, "$PWD/output/values.ini", sampler="test", outdir="test_run")

    dzp_early,dzp_late,dzm_early,dzm_late=[],[],[],[]

    if args.correct_dz:
        early = replace_dz_values("early",mean_values_dict)
        early_ref = replace_dz_values("none",mean_values_dict)
        export_values(early)
        dpath = "/share/des/disc3/samuroff/mpp/cosmosis/sws/ias/datavecs/real/mcal-gaussian/dvec-mcal-early.fits"
        ct.run_cosmosis(ini, "$PWD/output/values.ini", extra_args="2pt_like.data_file=%s fits_nz.nz_file=%s"%(dpath,dpath), sampler="test", outdir="test_run1")
        export_values(early_ref)
        ct.run_cosmosis(ini, "$PWD/output/values.ini", extra_args="2pt_like.data_file=%s fits_nz.nz_file=%s"%(dpath,dpath), sampler="test", outdir="test_run2")
        theta, dzp_early, dzm_early = get_shift()
        late = replace_dz_values("late",mean_values_dict)
        late_ref = replace_dz_values("none",mean_values_dict)
        export_values(late)
        dpath = "/share/des/disc3/samuroff/mpp/cosmosis/sws/ias/datavecs/real/mcal-gaussian/dvec-mcal-late.fits"
        ct.run_cosmosis(ini, "$PWD/output/values.ini", extra_args="2pt_like.data_file=%s fits_nz.nz_file=%s"%(dpath,dpath), sampler="test", outdir="test_run1")
        export_values(late_ref)
        ct.run_cosmosis(ini, "$PWD/output/values.ini", extra_args="2pt_like.data_file=%s fits_nz.nz_file=%s"%(dpath,dpath), sampler="test", outdir="test_run2")
        theta, dzp_late, dzm_late = get_shift()


    cornerplot("test_run", 
        "/share/des/disc3/samuroff/mpp/cosmosis/sws/ias/datavecs/real/mcal-gaussian/dvec-mcal-early.fits",
        "/share/des/disc3/samuroff/mpp/cosmosis/sws/ias/datavecs/real/mcal-gaussian/dvec-mcal-late.fits", dz=[theta,dzp_early,dzp_late,dzm_early,dzm_late])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--chain',"-c", type=str, action='store')
    parser.add_argument('--fkeep',"-k", default=1.0, type=float, action='store')
    parser.add_argument('--blind',"-b", action='store_true')
    parser.add_argument('--noia', action='store_true')
    parser.add_argument('--correct_dz', action='store_true')

    args = parser.parse_args()

    main(args)
