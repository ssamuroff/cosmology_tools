import numpy as np
import treecorr
import twopoint
import fitsio as fio
import healpy as hp
from numpy.lib.recfunctions import append_fields, rename_fields
from .stage import PipelineStage, TWO_POINT_NAMES
import os

CORES_PER_TASK=64

global_measure_2_point = None

def task(ijkl):
    i,j,k,l=ijkl
    global_measure_2_point.call_treecorr(i,j,k,l)

class Measure2Point(PipelineStage):
    name = "2pt"
    inputs = {
        "weight"        : ("nofz", "weight.npy")          ,
        "nz_source"     : ("nofz", "nz_source_zbin.npy")  ,
        "nz_lens"       : ("nofz", "nz_lens_zbin.npy")    ,
        "randoms"       : ("nofz", "randoms_zbin.npy")    ,
        "gold_idx"      : ("nofz", "gold_idx.npy")        ,
        "lens_idx"      : ("nofz", "lens_idx.npy")        ,
        "ran_idx"       : ("nofz", "ran_idx.npy")         ,

    }
    outputs = {
        "xip"    : "{rank}_xip.txt",
        "xim"    : "{rank}_xim.txt",
        "gammat" : "{rank}_gammat.txt",
        "wtheta" : "{rank}_wtheta.txt",
        "any"    : "*{rank}*",
    }

    def __init__(self,param_file):
        """
        Initialise object and load catalogs.
        """
        super(Measure2Point,self).__init__(param_file)

        def load_bin_weight_files(file, objid, sample=1):

            filename = self.input_path(file)
            gt = self.samples[sample-1]
            if (len(self.samples)>1):
                filename = filename.replace("nofz", "nofz_%s"%gt)
                if sample>1:
                    shape = getattr(self, "shape%d"%sample)
                else:
                    shape = self.shape
            else:
                shape = self.shape
            array    = np.load(filename)

            if np.array_equal(array[:,0], objid):
                if (file == "nz_source")&self.params['has_sheared']:
                    return array[:,1],array[:,2],array[:,3],array[:,4],array[:,5]
                else:
                    return array[:,1]
            else:
                import pdb ; pdb.set_trace()
                binning = np.ones((objid.size,6)) + 100
                binning[(shape["flags"]==0)] = array
                
                print array[:,0],len(array[:,0]),objid,len(objid)
                if not np.array_equal(binning[:,0], objid): import pdb ; pdb.set_trace() 
                else: return binning[:,1]

        import glob
        if (self.params['2pt_only'] == 'shear-shear')|(self.params['lensfile'] == 'None'):
            names = TWO_POINT_NAMES[:2]
        for n,name in enumerate(names):
            files = glob.glob(self.output_path('any').format(rank=name))
            for file in files:
                try:
                    os.remove(file)
                except OSError:
                    print 'OSerror removing old 2pt file'

        # Load data
        self.get_samples()

        for i,sample in enumerate(self.samples):
            print i,sample
            self.load_cat(suffix=i,split=sample)
            if (i>0):
                gold = getattr(self,"gold%d"%(i+1))
            else:
                gold = self.gold 

            source_binning = load_bin_weight_files("nz_source",gold['objid'], sample=i+1)
            weight         = load_bin_weight_files("weight",gold['objid'], sample=i+1)
            setattr(self, "source_binning%d"%(i+1), source_binning)
            if i==0:
                setattr(self, "weight", weight)
            else:
                setattr(self, "weight%d"%(i+1), weight)

        global global_measure_2_point
        global_measure_2_point = self

    def get_samples(self):
        if 'colour_bins' in self.params.keys():
            self.samples = self.params['colour_bins'].split()
        else:
            self.samples = ['all']

    def run(self):
        #This is a parallel job

        if hasattr(self.params['zbins'], "__len__"):
            self.zbins=len(self.params['zbins'])-1
        else:
            self.zbins=self.params['zbins']

        nbin = self.zbins
        ncbin = len(self.samples) 
        print "Will compute +/- correlations in %d redshift bins and %d colour bins (%d correlation)"%(nbin,ncbin, nbin*(nbin+1) * ncbin * ncbin   )
        all_calcs = [(i,j,k,l) for i in xrange(nbin) for j in xrange(nbin) for k in xrange(ncbin) for l in xrange(ncbin)]
        calcs=[]
        for i,j,k,l in all_calcs:
            if ((k==l) and (i<=j)) or  (k!=l):
                calcs.append((i,j,k,l))

        self.theta  = []
        self.xi     = []
        self.xierr  = []
        self.calc   = []

        if self.comm:
            from .mpi_pool import MPIPool
            pool = MPIPool(self.comm)
            pool.map(task, calcs)
            pool.close()
        else:
            map(task, calcs)

    def get_im3shape_type(self,shape,gtype):
        import pdb ; pdb.set_trace()
        if gtype.lower()=='bulge':
            return (shape['is_bulge']==1)
        else:
            return (shape['is_bulge']==0)

    def get_type_split(self, sample, shape, pz, pz_1p, pz_1m, pz_2p, pz_2m):
        if sample=='':
            sample = 0
        else:
            sample = int(sample)-1
        gt = self.samples[sample].split("_")

        print "Cutting to galaxy type",gt

        # If we need to split into bright/faint subsets, then include this in the mask first
        if len(gt)>1:
            subpop = gt[0]
            galaxy_type = gt[1]

            r = 30 - 2.5 * np.log10(shape["flux_r"])
            rp1,rm1 = 30 - 2.5 * np.log10(shape["flux_r_1p"]), 30 - 2.5 * np.log10(shape["flux_r_1m"])
            rp2,rm2 = 30 - 2.5 * np.log10(shape["flux_r_2p"]), 30 - 2.5 * np.log10(shape["flux_r_2m"])
            median_r = np.median(r[shape["flags"]==0])
            if (subpop.lower()=="bright"):
                mag_mask = (r<median_r), (rp1<median_r), (rm1<median_r), (rp2<median_r), (rm2<median_r)
            elif (subpop.lower()=="faint"):
                mag_mask = (r>median_r), (rp1>median_r), (rm1>median_r), (rp2>median_r), (rm2>median_r)
        else:
            galaxy_type=gt[0]
            mag_mask = [np.ones(pz["t_bpz"].size).astype(bool)]*5

        # Select either early- or late-type galaxies
        if galaxy_type=='early':
            mask = [(pz['t_bpz']<1), (pz_1p['t_bpz']<1), (pz_1m['t_bpz']<1), (pz_2p['t_bpz']<1), (pz_2m['t_bpz']<1)]
        elif galaxy_type=='late':
            mask =  [(pz['t_bpz']>=1), (pz_1p['t_bpz']>=1), (pz_1m['t_bpz']>=1), (pz_2p['t_bpz']>=1), (pz_2m['t_bpz']>=1)]
        elif galaxy_type=='all':
            mask = [np.ones(pz["t_bpz"].size).astype(bool)]*5
        
        final_mask = []
        for (m1,m2) in zip(mask,mag_mask):
            final_mask.append(m1&m2)

        return tuple(final_mask)

    def load_array(self,d,file):

        import time
        t0=time.time()

        if self.params['has_sheared'] & (file=='shapefile'):
            d['flags_1p'] = 'flags_select_1p'
            d['flags_1m'] = 'flags_select_1m'
            d['flags_2p'] = 'flags_select_2p'
            d['flags_2m'] = 'flags_select_2m'

        if self.params['pdf_type']=='pdf':
            keys = [key for key in d.keys() if (d[key] is not None)&(key is not 'pzstack')]
        else:
            keys = [key for key in d.keys() if (d[key] is not None)]

        if 'objid' in keys:
            dtypes = [('objid','i8')]
        else:
            raise ValueError('missing object id in '+file)
        dtypes += [(key,'f8') for key in keys if (key is not 'objid')]
        if self.params['pdf_type']=='pdf':
            dtypes  += [('pzstack_'+str(i),'f8') for i in range(len(self.params['pdf_z']))]

        fits = fio.FITS(self.params[file])[-1]
        array = fits.read(columns=[d[key] for key in keys])

        array = rename_fields(array,{v: k for k, v in d.iteritems()})

        if ('weight' not in array.dtype.names) & (file=='shapefile'):
            array = append_fields(array, 'weight', np.ones(len(array)), usemask=False)

        if self.params['pdf_type']=='pdf':
            for i in range(len(self.params['pdf_z'])):
                array['pzstack'+str(i)]=fits.read(columns=d['pzstack']+str(i))

        if np.any(np.diff(array['objid']) < 1):
            raise ValueError('misordered or duplicate ids in '+file) 

        return array

    def load_cat(self, split="all", suffix=0):

        if suffix==0:
            suffix = ""
        else:
            suffix = "%d"%(suffix+1)

        import time
        import importlib
        col = importlib.import_module('.'+self.params['dict_file'],'pipeline')

        t0 = time.time()

        gold      = self.load_array(col.gold_dict, 'goldfile')
        print 'Done goldfile',time.time()-t0, gold.dtype.names

        shape     = self.load_array(col.shape_dict, 'shapefile')

        print 'Done shapefile',time.time()-t0,shape.dtype.names
        pz        = self.load_array(col.pz_bin_dict, 'photozfile')
        print 'Done pzfile',time.time()-t0, pz.dtype.names
        if self.params['has_sheared']:
            pz_1p = self.load_array(col.pz_bin_dict, 'photozfile_1p')
            print 'Done pz1pfile',time.time()-t0, pz_1p.dtype.names
            pz_1m = self.load_array(col.pz_bin_dict, 'photozfile_1m')
            print 'Done pz1mfile',time.time()-t0, pz_1m.dtype.names
            pz_2p = self.load_array(col.pz_bin_dict, 'photozfile_2p')
            print 'Done pz2pfile',time.time()-t0, pz_2p.dtype.names
            pz_2m = self.load_array(col.pz_bin_dict, 'photozfile_2m')
            print 'Done pz2mfile',time.time()-t0, pz_2m.dtype.names
        pz_nofz   = self.load_array(col.pz_stack_dict, 'photozfile_nz')
        print 'Done pznofzfile',time.time()-t0, pz_nofz.dtype.names
        if self.params['lensfile'] != 'None':
            lens      = self.load_array(col.lens_dict, 'lensfile')
            print 'Done lensfile',time.time()-t0,lens.dtype.names
            lens_pz   = self.load_array(col.lens_pz_dict, 'lensfile')
            print 'Done lens_pzfile',time.time()-t0,lens_pz.dtype.names

        if 'm1' not in shape.dtype.names:
            shape = append_fields(shape, 'm1', shape['m2'], usemask=False)
        if 'm2' not in shape.dtype.names:
            shape = append_fields(shape, 'm2', shape['m1'], usemask=False)
        if self.params['oneplusm']==False:
            print 'converting m to 1+m'
            shape['m1'] = np.copy(shape['m1'])+1.
            shape['m2'] = np.copy(shape['m2'])+1.
        if 'c1' in shape.dtype.names:
            shape['e1'] -= shape['c1']
            shape['e2'] -= shape['c2']
            shape['c1'] = None
            shape['c2'] = None
        if self.params['flip_e2']==True:
            print 'flipping e2'
            shape['e2']*=-1
        if self.params['lensfile'] != 'None':
            if 'pzbin' not in lens_pz.dtype.names:
                lens_pz = append_fields(lens_pz, 'pzbin', lens_pz['pzstack'], usemask=False)
            if 'pzstack' not in lens_pz.dtype.names:
                lens_pz = append_fields(lens_pz, 'pzstack', lens_pz['pzbin'], usemask=False)

        if not ((len(gold)==len(shape))
            & (len(gold)==len(pz))
            & (len(gold)==len(pz_nofz))):
            raise ValueError('shape, gold, or photoz length mismatch')
        if self.params['has_sheared']:
            if not ((len(gold)==len(pz_1p))
                & (len(gold)==len(pz_1m))
                & (len(gold)==len(pz_2p))
                & (len(gold)==len(pz_2m))):
                raise ValueError('shape, gold, or photoz length mismatch')        
        if self.params['lensfile'] != 'None':
            if (len(lens)!=len(lens_pz)):
                raise ValueError('lens and lens_pz length mismatch') 

        if self.params['lensfile'] != 'None':
            keys = [key for key in col.ran_dict.keys() if (col.ran_dict[key] is not None)]
            fits = fio.FITS(self.params['randomfile'])[-1]

            dtypes=[(key,'f8') for key in keys]
            randoms =  np.empty(fits.read_header()['NAXIS2'], dtype = dtypes)
            for key in keys:
                randoms[key]=fits.read(columns=[col.ran_dict[key]])

        if self.params['test_run']==True:
            idx = np.load(self.input_path("gold_idx"))
            gol    = gold[idx]
            shape   = shape[idx]
            pz      = pz[idx]
            pz_nofz = pz_nofz[idx]
            if self.params['has_sheared']:
                pz_1p   = pz_1p[idx]
                pz_1m   = pz_1m[idx]
                pz_2p   = pz_2p[idx]
                pz_2m   = pz_2m[idx]
            if self.params['lensfile'] != 'None':
                idx = np.load(self.input_path("lens_idx"))
                lens    = lens[idx]
                lens_pz = lens_pz[idx]
                idx = np.load(self.input_path("ran_idx"))
                randoms = randoms[idx]

        if 'pzbin_col' in gold.dtype.names:
            mask = (gold['pzbin_col'] >= 0)
        else:
            mask = (pz['pzbin'] > self.params['zlims'][0]) & (pz['pzbin'] <= self.params['zlims'][1])
            if self.params['has_sheared']:
                mask_1p = (pz_1p['pzbin'] > self.params['zlims'][0]) & (pz_1p['pzbin'] <= self.params['zlims'][1])
                mask_1m = (pz_1m['pzbin'] > self.params['zlims'][0]) & (pz_1m['pzbin'] <= self.params['zlims'][1])
                mask_2p = (pz_2p['pzbin'] > self.params['zlims'][0]) & (pz_2p['pzbin'] <= self.params['zlims'][1])
                mask_2m = (pz_2m['pzbin'] > self.params['zlims'][0]) & (pz_2m['pzbin'] <= self.params['zlims'][1])

        if 'i3s_type' in self.params.keys():
            type_mask = self.get_im3shape_type(shape, self.params['i3s_type'])
            mask = mask & type_mask

        if (split!='all'):
            type_mask, tmask1p, tmask1m, tmask2p, tmask2m  = self.get_type_split(suffix, shape, pz, pz_1p, pz_1m, pz_2p, pz_2m )
            mask = mask & type_mask
            mask_1p = mask & tmask1p
            mask_1m = mask & tmask1m
            mask_2p = mask & tmask2p
            mask_2m = mask & tmask2m

        if 'flags' in shape.dtype.names:
            mask = mask & (shape['flags']==0)
            if self.params['has_sheared']:
                mask_1p = mask_1p & (shape['flags_1p']==0)
                mask_1m = mask_1m & (shape['flags_1m']==0)
                mask_2p = mask_2p & (shape['flags_2p']==0)
                mask_2m = mask_2m & (shape['flags_2m']==0)
        if 'flags' in pz.dtype.names:
            mask = mask & (pz['flags']==0)
            if self.params['has_sheared']:
                mask_1p = mask_1p & (pz_1p['flags']==0)
                mask_1m = mask_1m & (pz_1m['flags']==0)
                mask_2p = mask_2p & (pz_2p['flags']==0)
                mask_2m = mask_2m & (pz_2m['flags']==0)

        print 'hardcoded spt region cut'
        mask = mask & (shape['dec']<-35)
        if self.params['has_sheared']:
            mask_1p = mask_1p & (shape['dec']<-35)
            mask_1m = mask_1m & (shape['dec']<-35)
            mask_2p = mask_2p & (shape['dec']<-35)
            mask_2m = mask_2m & (shape['dec']<-35)

        if 'footprintfile' in self.params.keys():
            print 'cutting catalog to footprintfile'
            footmask = np.in1d(hp.ang2pix(4096, np.pi/2.-np.radians(shape['dec']),np.radians(shape['ra']), nest=False),
                               fio.FITS(self.params['footprintfile'])[-1].read()['HPIX'],assume_unique=False)
            mask = mask & footmask
            if self.params['has_sheared']:
                mask_1p = mask_1p & footmask
                mask_1m = mask_1m & footmask
                mask_2p = mask_2p & footmask
                mask_2m = mask_2m & footmask

        if self.params['has_sheared']:
            full_mask     = mask & mask_1p & mask_1m & mask_2p & mask_2m
            setattr(self, "mask%s"%suffix, mask[full_mask])
            setattr(self, "mask_1p%s"%suffix, mask_1p[full_mask])
            setattr(self, "mask_1m%s"%suffix, mask_1m[full_mask])
            setattr(self, "mask_2p%s"%suffix, mask_2p[full_mask])
            setattr(self, "mask_2m%s"%suffix, mask_2m[full_mask])
        else:
            full_mask  = mask
            setattr(self, "mask%s"%suffix, mask[full_mask])
        setattr(self, "gold%s"%suffix, gold[full_mask])
        setattr(self, "shape%s"%suffix, shape[full_mask])
        if self.params['lensfile']!='None': setattr(self, "randoms%s"%suffix,randoms)

        setattr(self,"pz%s"%suffix, None)
        setattr(self,"pz_1p%s"%suffix, None)
        setattr(self,"pz_1m%s"%suffix, None)
        setattr(self,"pz_2p%s"%suffix, None)
        setattr(self,"pz_2m%s"%suffix, None)
        setattr(self,"pz_nofz%s"%suffix, None)

        return



    def load_second_shapecat(self, split="all", suffix=0):

        import time
        import importlib
        col = importlib.import_module('.'+self.params['dict_file'],'pipeline')

        t0 = time.time()

        self.gold2      = self.load_array(col.gold_dict, 'goldfile')
        print 'Done goldfile',time.time()-t0,self.gold2.dtype.names

        self.shape2     = self.load_array(col.shape_dict, 'shapefile')
        if ("external_mask" in self.params.keys()):
            print "Loading additional selection mask from %s"%self.params['external_mask']
            select = fio.FITS(self.params['external_mask'])[-1].read()
            to_cut = np.invert(select["flag"].astype(bool))
            print "Cut removes %d/%d"%(self.shape2["flags"][to_cut].size, self.shape2["flags"].size)
            print "    applying to foreground galaxies"
            self.shape2["flags"][to_cut] = 16
            if "flags_1p" in self.shape2.dtype.names: self.shape2["flags_1p"][to_cut] = 16
            if "flags_1m" in self.shape2.dtype.names: self.shape2["flags_1m"][to_cut] = 16
            if "flags_2p" in self.shape2.dtype.names: self.shape2["flags_2p"][to_cut] = 16
            if "flags_2m" in self.shape2.dtype.names: self.shape2["flags_2m"][to_cut] = 16
            self.shape2["weight"][to_cut] = 0.0
        print 'Done shapefile',time.time()-t0,self.shape2.dtype.names
        self.pz_2        = self.load_array(col.pz_bin_dict, 'photozfile')
        print 'Done pzfile',time.time()-t0,self.pz_2.dtype.names
        if self.params['has_sheared']:
            self.pz_1p_2 = self.load_array(col.pz_bin_dict, 'photozfile_1p')
            print 'Done pz1pfile',time.time()-t0,self.pz_1p_2.dtype.names
            self.pz_1m_2 = self.load_array(col.pz_bin_dict, 'photozfile_1m')
            print 'Done pz1mfile',time.time()-t0,self.pz_1m_2.dtype.names
            self.pz_2p_2 = self.load_array(col.pz_bin_dict, 'photozfile_2p')
            print 'Done pz2pfile',time.time()-t0,self.pz_2p_2.dtype.names
            self.pz_2m_2 = self.load_array(col.pz_bin_dict, 'photozfile_2m')
            print 'Done pz2mfile',time.time()-t0,self.pz_2m_2.dtype.names
        self.pz_nofz_2   = self.load_array(col.pz_stack_dict, 'photozfile_nz')
        print 'Done pznofzfile',time.time()-t0,self.pz_nofz_2.dtype.names
        if self.params['lensfile'] != 'None':
            self.lens_2      = self.load_array(col.lens_dict, 'lensfile')
            print 'Done lensfile',time.time()-t0,self.lens_2.dtype.names
            self.lens_pz_2   = self.load_array(col.lens_pz_dict, 'lensfile')
            print 'Done lens_pzfile',time.time()-t0,self.lens_pz_2.dtype.names

        if 'm1' not in self.shape2.dtype.names:
            self.shape2 = append_fields(self.shape2, 'm1', self.shape2['m2'], usemask=False)
        if 'm2' not in self.shape.dtype.names:
            self.shape2 = append_fields(self.shape2, 'm2', self.shape2['m1'], usemask=False)
        if self.params['oneplusm']==False:
            print 'converting m to 1+m'
            self.shape2['m1'] = np.copy(self.shape2['m1'])+1.
            self.shape2['m2'] = np.copy(self.shape2['m2'])+1.
        if 'c1' in self.shape.dtype.names:
            self.shape2['e1'] -= self.shape2['c1']
            self.shape2['e2'] -= self.shape2['c2']
            self.shape2['c1'] = None
            self.shape2['c2'] = None
        if self.params['flip_e2']==True:
            print 'flipping e2'
            self.shape2['e2']*=-1
        if 'pzbin' not in self.lens_pz_2.dtype.names:
            self.lens_pz_2 = append_fields(self.lens_pz_2, 'pzbin', self.lens_pz_2['pzstack'], usemask=False)
        if 'pzstack' not in self.lens_pz_2.dtype.names:
            self.lens_pz_2 = append_fields(self.lens_pz_2, 'pzstack', self.lens_pz_2['pzbin'], usemask=False)

        if not ((len(self.gold2)==len(self.shape2))
            & (len(self.gold2)==len(self.pz_2))
            & (len(self.gold2)==len(self.pz_nofz_2))):
            raise ValueError('shape, gold, or photoz length mismatch')
        if self.params['has_sheared']:
            if not ((len(self.gold2)==len(self.pz_1p_2))
                & (len(self.gold2)==len(self.pz_1m_2))
                & (len(self.gold2)==len(self.pz_2p_2))
                & (len(self.gold2)==len(self.pz_2m_2))):
                raise ValueError('shape, gold, or photoz length mismatch')        
        if self.params['lensfile'] != 'None':
            if (len(self.lens_2)!=len(self.lens_pz_2)):
                raise ValueError('lens and lens_pz length mismatch') 

        if self.params['lensfile'] != 'None':
            keys = [key for key in col.ran_dict.keys() if (col.ran_dict[key] is not None)]
            fits = fio.FITS(self.params['randomfile'])[-1]

            dtypes=[(key,'f8') for key in keys]
            self.randoms_2 = np.empty(fits.read_header()['NAXIS2'], dtype = dtypes)
            for key in keys:
                self.randoms_2[key]=fits.read(columns=[col.ran_dict[key]])

        if self.params['test_run']==True:
            idx = np.load(self.input_path("gold_idx"))
            self.gold2    = self.gold2[idx]
            self.shape2   = self.shape2[idx]
            self.pz_2      = self.pz[idx]
            self.pz_nofz_2 = self.pz_nofz[idx]
            if self.params['has_sheared']:
                self.pz_1p_2   = self.pz_1p_2[idx]
                self.pz_1m_2   = self.pz_1m_2[idx]
                self.pz_2p_2   = self.pz_2p_2[idx]
                self.pz_2m_2   = self.pz_2m_2[idx]
            if self.params['lensfile'] != 'None':
                idx = np.load(self.input_path("lens_idx"))
                self.lens_2    = self.lens_2[idx]
                self.lens_pz_2 = self.lens_pz_2[idx]
                idx = np.load(self.input_path("ran_idx"))
                self.randoms_2 = self.randoms_2[idx]

        if 'pzbin_col' in self.gold2.dtype.names:
            mask = (self.gold2['pzbin_col'] >= 0)
        else:
            mask = (self.pz_2['pzbin'] > self.params['zlims'][0]) & (self.pz_2['pzbin'] <= self.params['zlims'][1])
            if self.params['has_sheared']:
                mask_1p = (self.pz_1p_2['pzbin'] > self.params['zlims'][0]) & (self.pz_1p_2['pzbin'] <= self.params['zlims'][1])
                mask_1m = (self.pz_1m_2['pzbin'] > self.params['zlims'][0]) & (self.pz_1m_2['pzbin'] <= self.params['zlims'][1])
                mask_2p = (self.pz_2p_2['pzbin'] > self.params['zlims'][0]) & (self.pz_2p_2['pzbin'] <= self.params['zlims'][1])
                mask_2m = (self.pz_2m_2['pzbin'] > self.params['zlims'][0]) & (self.pz_2m_2['pzbin'] <= self.params['zlims'][1])

        if (split!='all'):
            type_mask, tmask1p, tmask1m, tmask2p, tmask2m  = self.get_type_split(sample=int(suffix)-1)
            mask = mask & type_mask
            mask_1p = mask & tmask1p
            mask_1m = mask & tmask1m
            mask_2p = mask & tmask2p
            mask_2m = mask & tmask2m

        if 'flags' in self.shape2.dtype.names:
            mask = mask & (self.shape2['flags']==0)
            if self.params['has_sheared']:
                mask_1p = mask_1p & (self.shape2['flags_1p']==0)
                mask_1m = mask_1m & (self.shape2['flags_1m']==0)
                mask_2p = mask_2p & (self.shape2['flags_2p']==0)
                mask_2m = mask_2m & (self.shape2['flags_2m']==0)
        if 'flags' in self.pz_2.dtype.names:
            mask = mask & (self.pz_2['flags']==0)
            if self.params['has_sheared']:
                mask_1p = mask_1p & (self.pz_1p_2['flags']==0)
                mask_1m = mask_1m & (self.pz_1m_2['flags']==0)
                mask_2p = mask_2p & (self.pz_2p_2['flags']==0)
                mask_2m = mask_2m & (self.pz_2m_2['flags']==0)

        print 'hardcoded spt region cut'
        mask = mask & (self.shape2['dec']<-35)
        if self.params['has_sheared']:
            mask_1p = mask_1p & (self.shape2['dec']<-35)
            mask_1m = mask_1m & (self.shape2['dec']<-35)
            mask_2p = mask_2p & (self.shape2['dec']<-35)
            mask_2m = mask_2m & (self.shape2['dec']<-35)

        if 'footprintfile' in self.params.keys():
            print 'cutting catalog to footprintfile'
            footmask = np.in1d(hp.ang2pix(4096, np.pi/2.-np.radians(self.shape2['dec']),np.radians(self.shape2['ra']), nest=False),
                               fio.FITS(self.params['footprintfile'])[-1].read()['HPIX'],assume_unique=False)
            mask = mask & footmask
            if self.params['has_sheared']:
                mask_1p = mask_1p & footmask
                mask_1m = mask_1m & footmask
                mask_2p = mask_2p & footmask
                mask_2m = mask_2m & footmask

        if self.params['has_sheared']:
            full_mask     = mask & mask_1p & mask_1m & mask_2p & mask_2m
            self.mask_2     = mask[full_mask]
            self.mask_1p_2  = mask_1p[full_mask]
            self.mask_1m_2  = mask_1m[full_mask]
            self.mask_2p_2  = mask_2p[full_mask]
            self.mask_2m_2  = mask_2m[full_mask]
        else:
            full_mask  = mask
            self.mask_2  = mask[full_mask]
        self.gold2      = self.gold2[full_mask]
        self.shape2     = self.shape2[full_mask]

        self.pz_2 = None
        self.pz_1p_2 = None
        self.pz_1m_2 = None
        self.pz_2p_2 = None
        self.pz_2m_2 = None
        self.pz_nofz_2 = None

        return

    def load_data(self):

        import time
        import importlib
        col = importlib.import_module('.'+self.params['dict_file'],'pipeline')

        # def lowercase_array(arr):
        #     old_names = arr.dtype.names
        #     new_names = [name.lower() for name in old_names]
        #     renames   = dict(zip(old_names, new_names))
        #     return rename_fields(arr, renames)
        t0 = time.time()

        self.gold      = self.load_array(col.gold_dict, 'goldfile')
        print 'Done goldfile',time.time()-t0,self.gold.dtype.names
        self.shape     = self.load_array(col.shape_dict, 'shapefile')
        print 'Done shapefile',time.time()-t0,self.shape.dtype.names
        self.pz        = self.load_array(col.pz_bin_dict, 'photozfile')
        print 'Done pzfile',time.time()-t0,self.pz.dtype.names
        if self.params['has_sheared']:
            self.pz_1p = self.load_array(col.pz_bin_dict, 'photozfile_1p')
            print 'Done pz1pfile',time.time()-t0,self.pz_1p.dtype.names
            self.pz_1m = self.load_array(col.pz_bin_dict, 'photozfile_1m')
            print 'Done pz1mfile',time.time()-t0,self.pz_1m.dtype.names
            self.pz_2p = self.load_array(col.pz_bin_dict, 'photozfile_2p')
            print 'Done pz2pfile',time.time()-t0,self.pz_2p.dtype.names
            self.pz_2m = self.load_array(col.pz_bin_dict, 'photozfile_2m')
            print 'Done pz2mfile',time.time()-t0,self.pz_2m.dtype.names
        self.pz_nofz   = self.load_array(col.pz_stack_dict, 'photozfile_nz')
        print 'Done pznofzfile',time.time()-t0,self.pz_nofz.dtype.names
        if self.params['lensfile'] != 'None':
            self.lens      = self.load_array(col.lens_dict, 'lensfile')
            print 'Done lensfile',time.time()-t0,self.lens.dtype.names
            self.lens_pz   = self.load_array(col.lens_pz_dict, 'lensfile')
            print 'Done lens_pzfile',time.time()-t0,self.lens_pz.dtype.names

        if 'm1' not in self.shape.dtype.names:
            self.shape = append_fields(self.shape, 'm1', self.shape['m2'], usemask=False)
        if 'm2' not in self.shape.dtype.names:
            self.shape = append_fields(self.shape, 'm2', self.shape['m1'], usemask=False)
        if self.params['oneplusm']==False:
            print 'converting m to 1+m'
            self.shape['m1'] = np.copy(self.shape['m1'])+1.
            self.shape['m2'] = np.copy(self.shape['m2'])+1.
        if 'c1' in self.shape.dtype.names:
            self.shape['e1'] -= self.shape['c1']
            self.shape['e2'] -= self.shape['c2']
            self.shape['c1'] = None
            self.shape['c2'] = None
        if self.params['flip_e2']==True:
            print 'flipping e2'
            self.shape['e2']*=-1
        if self.params['lensfile']!='None':
            if 'pzbin' not in self.lens_pz.dtype.names:
                self.lens_pz = append_fields(self.lens_pz, 'pzbin', self.lens_pz['pzstack'], usemask=False)
            if 'pzstack' not in self.lens_pz.dtype.names:
                self.lens_pz = append_fields(self.lens_pz, 'pzstack', self.lens_pz['pzbin'], usemask=False)

        if not ((len(self.gold)==len(self.shape))
            & (len(self.gold)==len(self.pz))
            & (len(self.gold)==len(self.pz_nofz))):
            raise ValueError('shape, gold, or photoz length mismatch')
        if self.params['has_sheared']:
            if not ((len(self.gold)==len(self.pz_1p))
                & (len(self.gold)==len(self.pz_1m))
                & (len(self.gold)==len(self.pz_2p))
                & (len(self.gold)==len(self.pz_2m))):
                raise ValueError('shape, gold, or photoz length mismatch')        
        if self.params['lensfile'] != 'None':
            if (len(self.lens)!=len(self.lens_pz)):
                raise ValueError('lens and lens_pz length mismatch') 

        if self.params['lensfile'] != 'None':
            keys = [key for key in col.ran_dict.keys() if (col.ran_dict[key] is not None)]
            fits = fio.FITS(self.params['randomfile'])[-1]

            dtypes=[(key,'f8') for key in keys]
            self.randoms = np.empty(fits.read_header()['NAXIS2'], dtype = dtypes)
            for key in keys:
                self.randoms[key]=fits.read(columns=[col.ran_dict[key]])

        if self.params['test_run']==True:
            idx = np.load(self.input_path("gold_idx"))
            self.gold    = self.gold[idx]
            self.shape   = self.shape[idx]
            self.pz      = self.pz[idx]
            self.pz_nofz = self.pz_nofz[idx]
            if self.params['has_sheared']:
                self.pz_1p   = self.pz_1p[idx]
                self.pz_1m   = self.pz_1m[idx]
                self.pz_2p   = self.pz_2p[idx]
                self.pz_2m   = self.pz_2m[idx]
            if self.params['lensfile'] != 'None':
                idx = np.load(self.input_path("lens_idx"))
                self.lens    = self.lens[idx]
                self.lens_pz = self.lens_pz[idx]
                idx = np.load(self.input_path("ran_idx"))
                self.randoms = self.randoms[idx]

        if 'pzbin_col' in self.gold.dtype.names:
            mask = (self.gold['pzbin_col'] >= 0)
        else:
            mask = (self.pz['pzbin'] > self.params['zlims'][0]) & (self.pz['pzbin'] <= self.params['zlims'][1])
            if self.params['has_sheared']:
                mask_1p = (self.pz_1p['pzbin'] > self.params['zlims'][0]) & (self.pz_1p['pzbin'] <= self.params['zlims'][1])
                mask_1m = (self.pz_1m['pzbin'] > self.params['zlims'][0]) & (self.pz_1m['pzbin'] <= self.params['zlims'][1])
                mask_2p = (self.pz_2p['pzbin'] > self.params['zlims'][0]) & (self.pz_2p['pzbin'] <= self.params['zlims'][1])
                mask_2m = (self.pz_2m['pzbin'] > self.params['zlims'][0]) & (self.pz_2m['pzbin'] <= self.params['zlims'][1])

        if 'i3s_type' in self.params.keys():
            type_mask = self.get_im3shape_type(shape, self.params['i3s_type'])
            mask = mask & type_mask


        if ('type_split' in self.params.keys()):
            type_mask, tmask1p, tmask1m, tmask2p, tmask2m  = self.get_type_split(1)
            mask = mask & type_mask
            mask_1p = mask & tmask1p
            mask_1m = mask & tmask1m
            mask_2p = mask & tmask2p
            mask_2m = mask & tmask2m

        if 'flags' in self.shape.dtype.names:
            mask = mask & (self.shape['flags']==0)
            if self.params['has_sheared']:
                mask_1p = mask_1p & (self.shape['flags_1p']==0)
                mask_1m = mask_1m & (self.shape['flags_1m']==0)
                mask_2p = mask_2p & (self.shape['flags_2p']==0)
                mask_2m = mask_2m & (self.shape['flags_2m']==0)
        if 'flags' in self.pz.dtype.names:
            mask = mask & (self.pz['flags']==0)
            if self.params['has_sheared']:
                mask_1p = mask_1p & (self.pz_1p['flags']==0)
                mask_1m = mask_1m & (self.pz_1m['flags']==0)
                mask_2p = mask_2p & (self.pz_2p['flags']==0)
                mask_2m = mask_2m & (self.pz_2m['flags']==0)

        print 'hardcoded spt region cut'
        mask = mask & (self.shape['dec']<-35)
        if self.params['has_sheared']:
            mask_1p = mask_1p & (self.shape['dec']<-35)
            mask_1m = mask_1m & (self.shape['dec']<-35)
            mask_2p = mask_2p & (self.shape['dec']<-35)
            mask_2m = mask_2m & (self.shape['dec']<-35)

        if 'footprintfile' in self.params.keys():
            print 'cutting catalog to footprintfile'
            footmask = np.in1d(hp.ang2pix(4096, np.pi/2.-np.radians(self.shape['dec']),np.radians(self.shape['ra']), nest=False),
                               fio.FITS(self.params['footprintfile'])[-1].read()['HPIX'],assume_unique=False)
            mask = mask & footmask
            if self.params['has_sheared']:
                mask_1p = mask_1p & footmask
                mask_1m = mask_1m & footmask
                mask_2p = mask_2p & footmask
                mask_2m = mask_2m & footmask

        if self.params['has_sheared']:
            full_mask     = mask & mask_1p & mask_1m & mask_2p & mask_2m
            self.mask     = mask[full_mask]
            self.mask_1p  = mask_1p[full_mask]
            self.mask_1m  = mask_1m[full_mask]
            self.mask_2p  = mask_2p[full_mask]
            self.mask_2m  = mask_2m[full_mask]
        else:
            full_mask  = mask
            self.mask  = mask[full_mask]
        self.gold      = self.gold[full_mask]
        self.shape     = self.shape[full_mask]

        self.pz = None
        self.pz_1p = None
        self.pz_1m = None
        self.pz_2p = None
        self.pz_2m = None
        self.pz_nofz = None

        return
    

    def call_treecorr(self,i,j,k,l):
        """
        This is a wrapper for interaction with treecorr.
        """
        print "Running 2pt analysis on redshift pair {},{}, colour pair {},{}".format(i, j, k, l)
        #
        
        verbose=0
        # Cori value
        num_threads=CORES_PER_TASK

        theta,xip,xim,xiperr,ximerr = self.calc_shear_shear(i,j,k,l,verbose,num_threads)

        self.theta.append(theta)
        self.xi.append([xip,xim])
        self.xierr.append([xiperr,ximerr])
        self.calc.append((i,j,k,l))

        return 0

    def get_m(self,i,sample=1):
        suffix = ""*(sample==1)+ "%d"%sample*(sample!=1)

        shapes = getattr(self, "shape%s"%suffix)
        source_binning = getattr(self, "source_binning%s"%sample)
        if self.params['has_sheared']:
            mk1p = getattr(self, "mask_1p%s"%suffix)
            mk1m = getattr(self, "mask_1m%s"%suffix)
            mk2p = getattr(self, "mask_2p%s"%suffix)
            mk2m = getattr(self, "mask_2m%s"%suffix)
        else:
            m = getattr(self, "mask%s"%suffix)
            mk1p = np.ones(m.size).astype(bool)
            mk2p = np.ones(m.size).astype(bool)
            mk1m = np.ones(m.size).astype(bool)
            mk2m = np.ones(m.size).astype(bool)
        mk =  getattr(self, "mask%s"%suffix)

        if self.params['has_sheared']:
            mask = (source_binning[0] == i)
            mask_1p = (source_binning[1] == i)
            mask_1m = (source_binning[2] == i)
            mask_2p = (source_binning[3] == i)
            mask_2m = (source_binning[4] == i)
            m1   = np.mean(shapes['m1'][mask&mk])
            m2   = np.mean(shapes['m2'][mask&mk])
            m1   += (np.mean(shapes['e1'][mask_1p&mk1p]) - np.mean(shapes['e1'][mask_1m&mk1m])) / (2.*self.params['dg'])
            m2   += (np.mean(shapes['e2'][mask_2p&mk2p]) - np.mean(shapes['e2'][mask_2m&mk2m])) / (2.*self.params['dg'])
            m1   = m1*np.ones(len(mask))
            m2   = m2*np.ones(len(mask))
        else:
            mask = (source_binning == i)
            m1 = shapes['m1']
            m2 = shapes['m2']

        return m1, m2, mask&mk

    def calc_shear_shear(self,i,j,k,l,verbose,num_threads):

    	m1j,m2j,maskj = self.get_m(i, sample=l+1)
    	weightj = getattr(self,"weight"+("2"*l))
        shapej = getattr(self,"shape"+("2"*l))

        m1i,m2i,maski = self.get_m(i, sample=k+1)
        weighti = getattr(self,"weight"+("2"*k))
        shapei = getattr(self,"shape"+("2"*k))

        if self.params['has_sheared']:
            cat_i = treecorr.Catalog(g1=shapei['e1'][maski]/m1i[maski], g2=shapei['e2'][maski]/m2i[maski], w=weighti[maski], ra=shapei['ra'][maski], dec=shapei['dec'][maski], ra_units='deg', dec_units='deg')
        else:
            cat_i = treecorr.Catalog(g1=shapei['e1'][maski], g2=shapei['e2'][maski], w=weighti[maski], ra=shapei['ra'][maski], dec=shapei['dec'][maski], ra_units='deg', dec_units='deg')
            biascat_i = treecorr.Catalog(k=np.sqrt(shapei['m1'][maski]*shapei['m2'][maski]), w=weighti[maski], ra=shapei['ra'][maski], dec=shapei['dec'][maski], ra_units='deg', dec_units='deg')

        if self.params['has_sheared']:
            cat_j = treecorr.Catalog(g1=shapej['e1'][maskj]/m1j[maskj], g2=shapej['e2'][maskj]/m2j[maskj], w=weightj[maskj], ra=shapej['ra'][maskj], dec=shapej['dec'][maskj], ra_units='deg', dec_units='deg')
        else:
            cat_j = treecorr.Catalog(g1=shapej['e1'][maskj], g2=shapej['e2'][maskj], w=weightj[maskj], ra=shapej['ra'][maskj], dec=shapej['dec'][maskj], ra_units='deg', dec_units='deg')
            biascat_j = treecorr.Catalog(k=np.sqrt(shapej['m1'][maskj]*shapej['m2'][maskj]), w=weightj[maskj], ra=shapej['ra'][maskj], dec=shapej['dec'][maskj], ra_units='deg', dec_units='deg')

        gg = treecorr.GGCorrelation(nbins=self.params['tbins'], min_sep=self.params['tbounds'][0], max_sep=self.params['tbounds'][1], sep_units='arcmin', bin_slop=self.params['slop'], verbose=verbose,num_threads=num_threads)
        gg.process(cat_i,cat_j)
        if self.params['has_sheared']:
            norm = 1.
        else:
            kk = treecorr.KKCorrelation(nbins=self.params['tbins'], min_sep=self.params['tbounds'][0], max_sep=self.params['tbounds'][1], sep_units='arcmin', bin_slop=self.params['slop'], verbose=verbose,num_threads=num_threads)
            kk.process(biascat_i,biascat_j)
            norm = kk.xi

        theta=np.exp(gg.meanlogr)
        xip = gg.xip/norm
        xim = gg.xim/norm
        xiperr = ximerr = np.sqrt(gg.varxi)/norm

        return theta, xip, xim, xiperr, ximerr

    def write(self):
        """
        Write data to files.
        """
        if self.comm is None:
            rank = 0
        else:
            rank = self.comm.Get_rank()

        for n,name in enumerate(TWO_POINT_NAMES):
            if (n<2)&(self.params['2pt_only'].lower() not in [None,'shear-shear','all']):
                continue
            
            filename = self.output_path(name).format(rank=rank)
            f = None 
            
            for (theta,xi_data,ijkl) in zip(self.theta, self.xi, self.calc):
                i,j,k,l = ijkl
                if xi_data[n] is None:
                    continue
                if f is None:
                    f = open(filename, 'w')
                for theta,xi_theta in zip(theta, xi_data[n]):
                    f.write("{} {} {} {} {} {}\n".format(theta, i, j, k, l, xi_theta))
            if f is not None:
                f.close()
