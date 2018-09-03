import desdb
import fitsio as fi
import pandas as pd
import numpy as np

query="""
SELECT coadd_objects_id, mean_z, z_mc, template 
FROM NSEVILLA.PHOTOZ_TPL_Y1MOF_1 
"""


#
#q = """
#SELECT coadd_objects_id, mean_z, z_mc, template 
#FROM NSEVILLA.PHOTOZ_TPL_Y1MOF_1 
#ORDER BY coadd_objects_id 
#OFFSET %d ROWS FETCH NEXT 1000000 ROWS ONLY;"""
#
conn = desdb.Connection()
#base=1000000
#
#results=np.zeros(base*10, dtype=[("coadd_objects_id", int), ("mean_z", float), ("template", float), ("z_mc", float) ])
#
#ifile=0
#for i in xrange(1000):
#    print "batch %d"%i
#    query = q%(i*base)
#    print query, 
#    dat = conn.quick(query)
#    print "(%d)"%len(dat)
#
#    dat = pd.DataFrame(dat)
#    for name in results.dtype.names: results[name][i*base:(i+1)*base] = np.array(dat[name])
#    if results["coadd_objects_id"][-1]!=0:
#    	print "Saving results"
#        out = fi.FITS("PHOTOZ_TPL_Y1MOF_1-batch%d.fits"%j, "rw")
#        out.write(out)
#        out.close()
#        results=np.zeros(base*10, dtype=[("coadd_objects_id", int), ("mean_z", float), ("template", float), ("z_mc", float) ])
#        j+=1
#    del(dat)
print "running query..."
print query
dat = conn.quick(query)
dat = pd.DataFrame(dat)

print "Saving results..."
results=np.zeros(len(dat), dtype=[("coadd_objects_id", int), ("mean_z", float), ("template", float), ("z_mc", float) ])
for name in results.dtype.names: results[name] = np.array(dat[name])
out=fi.FITS("/share/des/disc6/samuroff/PHOTOZ_TPL_Y1MOF_1.fits", "rw")
print "done"
out.write(results)
out.close()

