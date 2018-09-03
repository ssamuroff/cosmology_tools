import numpy as np
import scipy, glob, argparse, yaml
import scipy.interpolate as sint
import fitsio as fi
import pylab as plt
import matplotlib
import os
from tools.cosmosis import cosmosis_tools as ct
import tools.emcee as mc
from string import Template as Tm


lims={}
lims["+"]={}
lims["-"]={}
lims["+"][(1,1)] = [7.195005, 250.0]
lims["+"][(2,1)] = [7.195005, 250.0]
lims["+"][(3,1)] = [5.715196, 250.0]
lims["+"][(4,1)] = [5.715196, 250.0]
lims["+"][(2,2)] = [4.539741, 250.0]
lims["+"][(3,2)] = [4.539741, 250.0]
lims["+"][(4,2)] = [4.539741, 250.0]
lims["+"][(3,3)] = [3.606045, 250.0]
lims["+"][(4,3)] = [3.606045, 250.0]
lims["+"][(4,4)] = [3.606045, 250.0]
lims["+"][(1,2)] = [7.195005, 250.0]
lims["+"][(1,3)] = [5.715196, 250.0]
lims["+"][(1,4)] = [5.715196, 250.0]
lims["+"][(2,3)] = [4.539741, 250.0]
lims["+"][(2,4)] = [4.539741, 250.0]
lims["+"][(3,4)] = [3.606045, 250.0]
lims["-"][(1,1)] = [90.579750, 250.0]
lims["-"][(2,1)] = [71.950053, 250.0]
lims["-"][(3,1)] = [71.950053, 250.0]
lims["-"][(4,1)] = [71.950053, 250.0]
lims["-"][(2,2)] = [57.151958, 250.0]
lims["-"][(3,2)] = [57.151958, 250.0]
lims["-"][(4,2)] = [45.397414, 250.0]
lims["-"][(3,3)] = [45.397414, 250.0]
lims["-"][(3,4)] = [45.397414, 250.0]
lims["-"][(4,4)] = [36.060448, 250.0]
lims["-"][(1,2)] = [71.950053, 250.0]
lims["-"][(1,3)] = [71.950053, 250.0]
lims["-"][(1,4)] = [71.950053, 250.0]
lims["-"][(2,3)] = [57.151958, 250.0]
lims["-"][(2,4)] = [45.397414, 250.0]
lims["-"][(4,3)] = [45.397414, 250.0]


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
    

def cornerplot(theory1, theory2, data1, data2, show_cuts=False):

    plt.switch_backend("pdf")
    plt.style.use("y1a1")
    matplotlib.rcParams["ytick.minor.visible"]=False
    matplotlib.rcParams["ytick.minor.width"]=0.1
    ni,nj=np.genfromtxt("%s/shear_xi/values.txt"%theory1[0]).T[2]
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
            xip_a,xim_a = get_real_spectra(i,j,data1)
            xip_b,xim_b = get_real_spectra(i,j,data2)

            posp = positions[(i+1,j+1,"+")]
            ax = plt.subplot(rows,cols,posp)
            ax.annotate("(%d, %d)"%(i+1,j+1), (3.8,5.8), textcoords='data', fontsize=11, )
            ax.yaxis.set_tick_params(which='minor', left='off', right='off')
            plt.ylim(-5,5)
            #plt.yscale("log")
            plt.xscale("log")
            if posp==19:
                plt.ylabel(r"$\Delta \xi_+(\theta)$", fontsize=12)
                plt.xlabel(r"$\theta$ / arcmin", fontsize=12)
                plt.yticks(fontsize=12)
            else:
                plt.yticks(visible=False)
            plt.xlim(2.2,270)
            plt.xticks([10,100],["10", "100"])
            #plt.yticks([-2,0,2,4,6,8],['-2', '0', '2', '4', '6', '8'])
            plt.axhline(0, color='k')


            if show_cuts:
                xlower,xupper = lims['+'][(i+1,j+1)]
                plt.axvspan(1e-6, xlower, color='gray',alpha=0.2)
                plt.axvspan(xupper, 500, color='gray',alpha=0.2)



            linestyles=['-',':','--','-']
            for iline,(t1, t2) in enumerate(zip(theory1,theory2)):
                
                xta,xip_theory_a,xim_theory_a = get_theory_spectra(i,j,t1)
                xtb,xip_theory_b,xim_theory_b = get_theory_spectra(i,j,t2)

                interp_xip_a = sint.interp1d(np.log10(xta), xip_theory_a, kind='cubic')
                interp_xip_b = sint.interp1d(np.log10(xtb), xip_theory_b, kind='cubic')

                xip_a_remapped = (interp_xip_a(np.log10(xip_a[0])))
                xip_b_remapped = (interp_xip_b(np.log10(xip_b[0])))
                erra = 1e4 * xip_a[0]*xip_a[1] * (xip_a[2]/xip_a[1])
                errb = 1e4 * xip_b[0]*xip_b[1] * (xip_b[2]/xip_b[1])


                plt.plot(xip_a[0], 1e5*(xip_a[1]-xip_a_remapped), ls=linestyles[iline], color="red")
                plt.plot(xip_b[0], 1e5*(xip_b[1]-xip_b_remapped), ls=linestyles[iline], color="royalblue")

            posm = positions[(i+1,j+1,"-")]
            ax = plt.subplot(rows,cols,posm)
            ax.annotate("(%d, %d)"%(i+1,j+1), (37,5.8), textcoords='data', fontsize=11, )
            ax.yaxis.set_tick_params(which='minor', left='off', right='off')
            plt.ylim(-5,5)
            if posm==30:
                ax.yaxis.set_label_position("right")
                ax.yaxis.set_ticks_position("right")
                plt.yticks(fontsize=12)
                plt.ylabel(r"$\Delta \xi_-(\theta)$", fontsize=12)
            else:
                plt.yticks(visible=False)
            #plt.yscale("log")
            plt.xscale("log")
            plt.xlim(2.2,270)
            plt.xticks([10,100],["10", "100"])
            #plt.yticks([-2,0,2,4,6,8],['-2', '0', '2', '4', '6', '8'])
            plt.axhline(0, color='k')

            if show_cuts:
                xlower,xupper = lims['-'][(i+1, j+1)]
                plt.axvspan(1e-6, xlower, color='gray',alpha=0.2)
                plt.axvspan(xupper, 500, color='gray',alpha=0.2)

            for iline,(t1, t2) in enumerate(zip(theory1,theory2)):
                
                xta,xip_theory_a,xim_theory_a = get_theory_spectra(i,j,t1)
                xtb,xip_theory_b,xim_theory_b = get_theory_spectra(i,j,t2)

                interp_xim_a = sint.interp1d(np.log10(xta), xim_theory_a, kind='cubic')
                interp_xim_b = sint.interp1d(np.log10(xtb), xim_theory_b, kind='cubic')

                xim_a_remapped = (interp_xim_a(np.log10(xim_a[0])))
                xim_b_remapped = (interp_xim_b(np.log10(xim_b[0])))
                erra = 1e4 * xim_a[0]*xim_a[1] * (xim_a[2]/xim_a[1])
                errb = 1e4 * xim_b[0]*xim_b[1] * (xim_b[2]/xim_b[1])


                plt.plot(xim_a[0], 1e5*(xim_a[1]-xim_a_remapped), ls=linestyles[iline], color="red")
                plt.plot(xim_b[0], 1e5*(xim_b[1]-xim_b_remapped), ls=linestyles[iline], color="royalblue")


    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig("/home/samuroff/tmp3.pdf")



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
        "bias_4" : 0.0},
    "lens_photoz_errors":
        {
        "bias_1" : 0.0,
        "bias_2" : 0.0,
        "bias_3" : 0.0,
        "bias_4" : 0.0},
    "bias_parameters":
        {
        "b_1" : 0.0,
        "b_2" : 0.0,
        "b_3" : 0.0,
        "b_4" : 0.0}}

    for i, name in enumerate(mean_vals["name"]):
        if (name=="post") or (name=="weight"):
            continue
        section, param = name.split("--")
        print section, param,
        value = mean_vals["mean"][i]
        params[section][param] = value
        if blind:
            print "XXX"
        else:
            print value

    return params

def export_values(mean_values_dict, export_to='output/values.ini'):
    vals = Tm(open("/home/samuroff/local/python/lib/python2.7/site-packages/tools/cosmosis/values_template").read())
    vals_all = {}
    for name1 in mean_values_dict.keys():
        for  name2 in mean_values_dict[name1].keys():
            vals_all[name2] = mean_values_dict[name1][name2]

    values_txt = vals.substitute(vals_all)

    outfile = open(export_to, "wa")
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
   
    cornerplot(args.lines1, args.lines2, 
        args.data1,
        args.data2, show_cuts=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--lines1',"-l1", nargs='+', action='store', default='none')
    parser.add_argument('--lines2',"-l2", nargs='+', action='store', default='none')
    parser.add_argument('--data1',"-d1", type=str, action='store')
    parser.add_argument('--data2',"-d2", type=str, action='store')

    args = parser.parse_args()


    main(args)
