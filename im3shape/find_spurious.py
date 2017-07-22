import glob
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
import sys

fsize=10   # width of features we are looking for in y pix 
csize=100  # width of range in y used to estimate / subtract background
dim=10000  # size of tiles in pix

def my_kde_left(y): # boxcar kernel density estimate with compensation of left
    ys=np.zeros(dim)
    y=np.sort(y)
    fc=0 # first and last index of compensation range
    lc=0
    ff=0 # first and last index of positive part of filter
    lf=0
    
    for i in range(dim):
        while(y[fc]<i-fsize/2-csize):
            fc += 1
        while(y[lc]<i-fsize/2):
            lc += 1
        ff=lc
        while(y[lf]<i+fsize/2 and lf<len(y)-1):
            lf += 1
        ys[i] = (lf-ff)-(lc-fc)*fsize/csize
    return ys
     
def my_kde_right(y): # boxcar kernel density estimate with compensation on right
    ys=np.zeros(dim)
    y=np.sort(y)
    fc=0 # first and last index of compensation range
    lc=0
    ff=0 # first and last index of positive part of filter
    lf=0
    
    for i in range(dim):
        while(y[fc]<i+fsize/2 and fc<len(y)-1):
            fc += 1
        while(y[lc]<i+fsize/2+csize and lc<len(y)-1):
            lc += 1
        lf=fc
        while(y[ff]<i-fsize/2):
            ff += 1
        ys[i] = (lf-ff)-(lc-fc)*fsize/csize
    return ys
        

if(len(sys.argv)!=2):
  print("syntax:",sys.argv[0],"[coadd re-run SExtractor catalog]")
  sys.exit(1)

f=fits.open(sys.argv[1])
print(sys.argv[1])

yl=my_kde_left(f[1].data["Y_IMAGE"])
yr=my_kde_right(f[1].data["Y_IMAGE"])
x=np.arange(0,dim)

#plt.plot(x,yl,c='b',alpha=0.5)
#plt.plot(x,yr,c='g',alpha=0.5)
#plt.savefig(sys.argv[1]+"_fhist.eps")

ybad=x[np.logical_and(yl>25,yr>25)] # threshold to cut; 20 gives a lot of false positives

plt.switch_backend("agg")
plt.scatter(f[1].data["X_IMAGE"],f[1].data["Y_IMAGE"],marker="o",s=0.02,alpha=0.5)
plt.scatter(np.zeros(len(ybad))-500,ybad,c='r')
plt.savefig(sys.argv[1]+"_fp.png",dpi=600)

tile=os.path.basename(sys.argv[1])[:12]

out=open("bad_rows", "a")
for r in ybad: out.write("%s %d \n"%(tile,r))

out.close()
# anything with Y_IMAGE coordinate in any of the entries of ybad +-fsize/2 is possibly spurious --
