import tools.diagnostics as di
import tools.plots as pl
import fitsio as fi
import pylab as plt
import os, yaml, argparse
import tools.shapes as s
import tools.nbc as cal
import tools.arrays as arr
from tools.im3shape import calibrate_all_tomographic as ct
from tools.im3shape import calibrate_all as ca
import copy
from tools.im3shape.calibration import calibrate_y1 as cy1
config = yaml.load(open("/home/samuroff/local/python/lib/python2.7/site-packages/tools/im3shape/calibration/config/fiducial-y1-unblinded.yaml"))
hoopoe, weights, y1v2_noblind = cy1.setup(True, True, config)

plt.style.use("y1a1")
plt.switch_backend("pdf")
fig=plt.figure(1)
pl.footprint_sub(y1v2_noblind.res["ra"],y1v2_noblind.res["dec"],10,10,4096,fig)
plt.savefig('/home/samuroff/shear_pipeline/plot_dump/y1a1-im3shape-footprint.pdf')
plt.close()
