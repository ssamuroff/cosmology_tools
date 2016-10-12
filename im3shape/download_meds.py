import numpy as np
import galsim
import fitsio
import yaml
import glob, argparse, os
import tools.shapes as s


def main():
    global args

    parser = argparse.ArgumentParser(add_help=True)
    parser.add_argument('-c', '--configuraion', type=str, default=".", help='hoopoe config file we want to run.')
    parser.add_argument('--meds', action="store_true", help='download MEDS files.')
    parser.add_argument('--images', action="store_true", help='download input image files for each tile.')

    args = parser.parse_args()

    config = yaml.load(open(args.configuraion))
    tilelist = yaml.load(open(config["coadd_tilelist"]))["tile_ids"]
    bands = config["bands"]

    try:
        local_path = config["meds_dir"]
    except:
        local_path = "%s/%s/final/"%(config["desdata"],config["meds_id"])

    print "Will put MEDS files in directory %s"%local_path


    for band in bands:
        nersc_path = "/project/projectdirs/des/wl/desdata/users/cogs/meds/%s/"%band
        for i, tile in enumerate(tilelist):

            filename = "%s*-%s-*.fits.fz"%(tile,band)
            print i, filename

            if args.meds:
                command = "rsync -urv sws@edison.nersc.gov:%s/%s %s/"%(nersc_path,filename,local_path)
                os.system(command)

            if args.images:
    
                meds_path = glob.glob("%s/%s"%(local_path,filename))[0]
                meds = fitsio.FITS(meds_path)
                image_table = meds["image_info"].read()

                for im in image_table["image_path"]:
                    path = os.path.dirname(im.strip().replace("/astro/u/astrodat/data/DES/", config["desdata"]))
                    remote_path = im.strip().replace("/astro/u/astrodat/data/DES/","/global/cscratch1/sd/sws/data/")


                    os.system("mkdir -p %s"%path)
                    if "coadd" in remote_path:
                        command = "scp -r sws@edison.nersc.gov:%s/%s*_%s* %s/"%(os.path.dirname(remote_path),tile,band,path)
                    else:
                        command = "scp -r sws@edison.nersc.gov:%s %s/"%(remote_path,path)
                    
                    os.system(command)


print "Done"


main()
