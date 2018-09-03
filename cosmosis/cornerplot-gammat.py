import numpy as np
import scipy, glob, argparse, yaml
import fitsio as fi
import pylab as plt
import matplotlib
import os
from tools.cosmosis import cosmosis_tools as ct
import tools.emcee as mc
from string import Template as Tm


lims={}

lims[(1,1)] = [64.0, 250.0]
lims[(2,1)] = [40.0, 250.0]
lims[(3,1)] = [30.0, 250.0]
lims[(4,1)] = [24.0, 250.0]
lims[(2,2)] = [40.0, 250.0]
lims[(3,2)] = [30.0, 250.0]
lims[(4,2)] = [24.0, 250.0]
lims[(3,3)] = [30.0, 250.0]
lims[(4,3)] = [24.0, 250.0]
lims[(4,4)] = [24.0, 250.0]
lims[(1,2)] = [64.0, 250.0]
lims[(1,3)] = [64.0, 250.0]
lims[(1,4)] = [64.0, 250.0]
lims[(2,3)] = [40.0, 250.0]
lims[(2,4)] = [40.0, 250.0]
lims[(3,4)] = [30.0, 250.0]

lims[(1,5)] = [5.715196, 250.0]
lims[(2,5)] = [7.195005, 250.0]
lims[(3,5)] = [7.195005, 250.0]
lims[(4,5)] = [7.195005, 250.0]

lims[(5,1)] = [21.0, 250.0]
lims[(5,2)] = [21.0, 250.0]
lims[(5,3)] = [21.0, 250.0]
lims[(5,4)] = [21.0, 250.0]


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
    theta = np.loadtxt("%s/galaxy_shear_xi/theta.txt"%filename)
    theta = theta * 60. * 180. / np.pi
    gammat = np.loadtxt("%s/galaxy_shear_xi/bin_%d_%d.txt"%(filename,j+1,i+1))
    return theta, gammat

def get_real_spectra(i,j,fits,error=True):
    if fits is None:
        return [],[],[]
    else:
        selectp = (fits["gammat"]["BIN1"][:]==j+1) & (fits["gammat"]["BIN2"][:]==i+1)

        xp = fits["gammat"]["ANG"][:][selectp]

        y = fits["gammat"]["VALUE"][:][selectp]

        if error:
            cov = fits["COVMAT"][:,:]
            start,end = fits["covmat"].read_header()["STRT_2"], fits["covmat"].read_header()["STRT_3"]
            dx = end - start

            nx = xp.size

            covp = cov[start:end,start:end]

            i0 = nx*(i + j)
            err = np.diag(covp[i0:(i0+nx),i0:(i0+nx)])
            err = np.sqrt(err)

        else:
            err = None

        return (xp,y,err)

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
    ni,nj=np.genfromtxt("%s/galaxy_shear_xi/values.txt"%theory1).T[2]
    ni = int(ni)
    nj = int(nj)
    if data1 is not None:
        data1 = fi.FITS(data1)
    if data2 is not None:
        data2 = fi.FITS(data2)

    fig = plt.figure(0)

    rows, cols = nj, ni

    count = 0


    for i in xrange(nj):
        for j in xrange(ni):
            count+=1

            print i,j
            xta, gammat_theory_a = get_theory_spectra(i,j,theory1)
            xtb, gammat_theory_b = get_theory_spectra(i,j,theory2)
            gammat_a  = get_real_spectra(i,j,data1)
            gammat_b  = get_real_spectra(i,j,data2)

            ax = plt.subplot(rows,cols, count)
            ax.annotate("(%d, %d)"%(i+1,j+1), (3.8,3), textcoords='data', fontsize=11, )
            ax.yaxis.set_tick_params(which='minor', left='off', right='off')
            plt.ylim(-2,4)
            #plt.yscale("log")
            plt.xscale("log")
            if (j==0):
                fig.text(0.07,0.5, r"$\theta \gamma_t(\theta)$ / $10^{-2}$ arcmin", ha='center', va='center', rotation='vertical', fontsize=18)
            if (i==nj-1):
                plt.xlabel(r"$\theta$ / arcmin", fontsize=12)
            plt.yticks(fontsize=12)

            if not (j==0):
                plt.yticks(visible=False)
            plt.xlim(2.2,270)
            plt.xticks([10,100],["10", "100"])
            plt.yticks([-1,0,1,2,3],['-1', '0', '1', '2', '3'])
            plt.axhline(0, color='k')


            if show_cuts:
                xlower,xupper = lims[(j+1,i+1)]
                plt.axvspan(1e-6, xlower, color='gray',alpha=0.2)
                plt.axvspan(xupper, 500, color='gray',alpha=0.2)

            #import pdb ; pdb.set_trace()
            erra = 1e2 * gammat_a[0]*gammat_a[1] * (gammat_a[2]/gammat_a[1])
            plt.errorbar(gammat_a[0], 1e2 * gammat_a[0]*gammat_a[1], yerr=erra, ls="none", marker="D", markersize=4,  ecolor="red", markeredgecolor="red", markerfacecolor="none")
            errb = 1e2 * gammat_b[0]*gammat_b[1] * (gammat_b[2]/gammat_b[1])
            plt.errorbar(gammat_b[0], 1e2 * gammat_b[0]*gammat_b[1], yerr=errb, ls="none", marker="o", markersize=4, ecolor="royalblue", markeredgecolor="royalblue", markerfacecolor="royalblue")

            plt.plot(xta, 1e2*xta*gammat_theory_a,color="k", ls='--')
            plt.plot(xtb, 1e2*xtb*gammat_theory_b,color="k")


    plt.subplots_adjust(hspace=0,wspace=0)
    plt.savefig("/home/samuroff/tmp4.pdf")



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
    vals = Tm(open("/home/samuroff/local/python/lib/python2.7/site-packages/tools/cosmosis/values_template-3x2pt").read())
    vals_all = {}
    for name1 in mean_values_dict.keys():
        for  name2 in mean_values_dict[name1].keys():
            vals_all[name2] = mean_values_dict[name1][name2]

    values_txt = vals.substitute(vals_all)

    outfile = open(export_to, "wa")
    outfile.write(values_txt)
    outfile.close()


def main(args):
   
    if (args.chain1!='none') and (not args.disable_cosmosis):
         # postprocess the specified chains
        os.system("mkdir -p output1/")
        print "Copying chain", args.chain1
        os.system("cp %s output1/chain1.txt"%args.chain1)
        ct.postprocess("output1/chain1.txt", "output1")

        # Read and convert the mean values to an ini file cosmosis can run with
        mean_values_dict = parse_values("output1", args, blind=args.blind)

        # Write the marginalised parameter means to an ini file 
        export_values(mean_values_dict, export_to="output1/values.ini")

        # Finally call cosmosis
        ini = "/home/samuroff/local/python/lib/python2.7/site-packages/tools/cosmosis/ini/des-test-3x2pt.ini"
        os.putenv("DATAFILE", '/share/des/disc3/samuroff/mpp/cosmosis/sws/ias/datavecs/real/mcal-nongaussian/dvec-multiprobe-mcal-early-NG.fits')
        ct.run_cosmosis(ini, "$PWD/output1/values.ini", sampler="test", outdir="test_run1")
        


    if (args.chain2!='none') and (not args.disable_cosmosis):
        os.system("mkdir -p output2/")
        print "Copying chain", args.chain2
        os.system("cp %s output2/chain2.txt"%args.chain2)
        ct.postprocess("output2/chain2.txt", "output2")

        mean_values_dict = parse_values("output2", args, blind=args.blind)

        export_values(mean_values_dict, export_to="output2/values.ini")

        ini = "/home/samuroff/local/python/lib/python2.7/site-packages/tools/cosmosis/ini/des-test-3x2pt.ini"
        os.putenv("DATAFILE", '/share/des/disc3/samuroff/mpp/cosmosis/sws/ias/datavecs/real/mcal-nongaussian/dvec-multiprobe-mcal-late-NG.fits')
        ct.run_cosmosis(ini, "$PWD/output2/values.ini", sampler="test", outdir="test_run2")
        

    cornerplot("test_run1", "test_run2", 
        args.data1,
        args.data2, show_cuts=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--chain1',"-c1", type=str, action='store', default='none')
    parser.add_argument('--chain2',"-c2", type=str, action='store', default='none')
    parser.add_argument('--data1',"-d1", type=str, action='store')
    parser.add_argument('--data2',"-d2", type=str, action='store')
    parser.add_argument('--blind',"-b", action='store_true')
    parser.add_argument('--disable_cosmosis', action='store_true')

    args = parser.parse_args()

    main(args)
