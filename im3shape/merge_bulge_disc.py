from astropy.table import Table, vstack, Column
# astropy FTW
import numpy as np
import glob, os
from astropy.io import fits as pyfits




def merge(bulge_main, bulge_epoch, disc_main, disc_epoch, bord_main, bord_epoch):
    disc_ids, bulge_ids = merge_main(bulge_main, disc_main, bord_main)
    merge_epoch(bulge_epoch, disc_epoch, bord_epoch, disc_ids, bulge_ids)

def merge_epoch(bulge_filename, disc_filename, bord_filename, disc_ids, bulge_ids):
    cat_b = Table.read(bulge_filename)
    cat_d = Table.read(disc_filename)
    w_b = np.in1d(cat_b['coadd_objects_id'], bulge_ids)
    w_d = np.in1d(cat_d['coadd_objects_id'], disc_ids)
    cat_b = cat_b[w_b]
    cat_d = cat_d[w_d]
    cat = vstack([cat_b, cat_d])
    cat.sort('coadd_objects_id')
    print "Writing merged epoch file {},  using {} from disc and {} from bulge".format(bord_filename, len(cat_d), len(cat_b))
    cat.write(bord_filename)



def merge_main(bulge_filename, disc_filename, bord_filename):
    if os.path.exists(bord_filename):
        print "Already done: ", bord_filename, len(bulge_files)
        return

    assert os.path.split(disc_filename)[1]==os.path.split(bulge_filename)[1]
    cat_b = Table.read(bulge_filename)
    cat_d = Table.read(disc_filename)
    print bulge_filename,len(cat_b)
    print disc_filename,len(cat_d)
    import pdb ; pdb.set_trace()

    if cat_b["covmat_0_1"][0].dtype is not float:
        nparam = cat_b["nparam_varied"]
        for i in xrange(nparam):
            for j in xrange(nparam):
                cat_b["covmat_%d_%d"%(i,j)]=cat_b["covmat_%d_%d"%(i,j)].astype(float)

    import pdb ; pdb.set_trace()
    
    
    # need to match coadd_object_ids here
    intersect = np.intersect1d(cat_b['coadd_objects_id'],cat_d['coadd_objects_id'])
    
    select_b = np.where(np.in1d(cat_b['coadd_objects_id'],intersect))[0]
    select_d = np.where(np.in1d(cat_d['coadd_objects_id'],intersect))[0]
    
    ndiscard_bulge =  len(cat_b) - len(select_b)
    ndiscard_disc =   len(cat_d) - len(select_d)
    cat_b = cat_b[select_b]
    cat_d = cat_d[select_d]
    print 'len intersect b=',len(cat_b),'len intersect d=',len(cat_d),'ndiscard_disc d=',ndiscard_disc, 'ndiscard_b=', ndiscard_bulge
    
    assert (cat_b['coadd_objects_id'] == cat_d['coadd_objects_id']).all()

    bulge_good = (cat_b['info_flag']==0) & (cat_b['error_flag']==0)
    disc_good = (cat_d['info_flag']==0) & (cat_d['error_flag']==0)
    # so I can understand the next bit
    bulge_bad = ~bulge_good
    disc_bad = ~disc_good

    bulge_better = (cat_b['likelihood']-cat_d['likelihood']) > 0
    disc_better = (cat_b['likelihood']-cat_d['likelihood']) < 0
    
    
    # SO in this bit we allow the info_flag to override
    # the likelihood in deciding which is better?
    # do we want that?
    
    select_b = (bulge_good&disc_bad) | (bulge_good&disc_good&bulge_better) | (bulge_bad&disc_bad&bulge_better)
    select_d = (disc_good&bulge_bad) | (bulge_good&disc_good&disc_better) | (bulge_bad&disc_bad&disc_better)
    
    cat_b = cat_b[select_b]
    cat_d = cat_d[select_d]

    disc_ids = cat_d['coadd_objects_id']
    bulge_ids = cat_b['coadd_objects_id']


    cat_final = vstack([cat_b, cat_d])
    cat_final.sort('coadd_objects_id')

    cat_final.add_column(Column(name='is_bulge', data=(cat_final['bulge_a']>0).astype(int)))
        
    n_total = len(cat_final)
    n_common = np.sum(bulge_good & disc_good)
    n_common_b = np.sum(bulge_good & disc_good & bulge_better)
    n_common_d = np.sum(bulge_good & disc_good & disc_better)
    n_only_b = np.sum(bulge_good & (~disc_good))
    n_only_d = np.sum(disc_good & (~bulge_good))  
    n_check = n_only_b + n_only_d +  n_common_b + n_common_d    
    
    cat_final.write(bord_filename)
    print '%s  n_total=%d n_common=%d n_common_b=%d n_common_d=%d n_only_b=%d n_only_d=%d n_check=%d' % (bord_filename,n_total,n_common,n_common_b,n_common_d,n_only_b,n_only_d,n_check)

    return disc_ids, bulge_ids


import sys
if '--mpi' in sys.argv:
    from mpi4py.MPI import COMM_WORLD as world
    rank = world.Get_rank()
    size = world.Get_size()
else:
    rank = 0
    size = 1

def main():

    bulge_path = 'bulge/main/DES*'
    disc_path = 'disc/main/DES*'

    print "Running merge script."

    bulge_files = glob.glob(bulge_path)
    disc_files = glob.glob(disc_path)

    bulge_tiles = [os.path.split(b)[1] for b in bulge_files]
    disc_tiles = [os.path.split(b)[1] for b in disc_files]

    print bulge_tiles[:10]
    print disc_tiles[:10]

    bulge_files = [b for i,b in enumerate(bulge_files) if bulge_tiles[i] in disc_tiles]
    disc_files = [b for i,b in enumerate(disc_files) if disc_tiles[i] in bulge_tiles]


    bulge_files.sort()
    disc_files.sort()

     
    ifiles = 0


    for i,(bulge_filename, disc_filename) in enumerate(zip(bulge_files, disc_files)):
        if i%size!=rank: continue
        try:
            bord_filename = bulge_filename.replace('bulge','bord')
            merge_main(bulge_filename, disc_filename, bord_filename)
        except Exception as error:
            print error
main()
