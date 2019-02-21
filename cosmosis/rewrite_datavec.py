import numpy as np
import fitsio as fi
import os
import argparse
import scipy.interpolate as interp

class output:
    def __init__(self, data):
        filename = data.replace('.fits','-mock.fits')
        print 'Generating new file: %s'%filename

        os.system('cp %s %s'%(data,filename))

        self.fits = fi.FITS(filename,'rw')

    def get_path(self, path, correlation):
        if correlation=='xip':
            pth = '%s/shear_xi/'%path + 'xiplus_%d_%d.txt'
            xpth = '%s/shear_xi/theta.txt'%path
        elif correlation=='xim':
            pth = '%s/shear_xi/'%path + 'ximinus_%d_%d.txt'
            xpth = '%s/shear_xi/theta.txt'%path
        elif correlation=='gammat':
            pth = '%s/galaxy_shear_xi/'%path + 'bin_%d_%d.txt'
            xpth = '%s/galaxy_shear_xi/theta.txt'%path
        elif correlation=='wtheta':
            pth = '%s/galaxy_xi/'%path + 'bin_%d_%d.txt'
            xpth = '%s/galaxy_xi/theta.txt'%path

        return xpth, pth

    def find_data(self, path, corrs=['xip','xim','gammat','wtheta']):
        
        for corr in corrs:
            print corr
            dvec = self.fits[corr].read()
            out = np.zeros(dvec['VALUE'].size, dtype=dvec.dtype)

            for name in dvec.dtype.names:
                out[name] = dvec[name]

            xpth, pth = self.get_path(path, corr)
            x = np.loadtxt(xpth)  * 180 / np.pi * 60
            done = []

            for j in dvec['BIN1']:
                for i in dvec['BIN2']:
                    if (i,j) in done:
                        continue
                    #print i, j
                    xdat = dvec['ANG'][(dvec['BIN1']==j) & (dvec['BIN2']==i)]
                    if os.path.exists(pth%(j,i)):
                        simulated = np.loadtxt(pth%(j,i))
                    elif ('xi' in corr) or ('wtheta' in corr) :
                        simulated = np.loadtxt(pth%(i,j))
                    else:
                        raise ValueError()

                    if (simulated.min()>0):
                        interpolator = interp.interp1d(np.log10(x), np.log10(simulated))
                        interpolated = 10**(interpolator(np.log10(xdat)))
                    else:
                        interpolator = interp.interp1d(np.log10(x), simulated)
                        interpolated = interpolator(np.log10(xdat))

                    out['VALUE'][(dvec['BIN1']==j) & (dvec['BIN2']==i)] = interpolated
                    done.append((i,j))

            print 'Saving to FITS file'
            self.fits[corr].write(out)
        self.fits.close()

        print 'Done'


def main(args):
    outfits = output(args.data)

    outfits.find_data(args.simulated, corrs=['xip','xim','gammat','wtheta'])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('--check',"-c", action='store_true', default=False)
    parser.add_argument('--simulated',"-s", action='store', default='none')
    parser.add_argument('--data',"-d", type=str, action='store')

    args = parser.parse_args()


    main(args)
