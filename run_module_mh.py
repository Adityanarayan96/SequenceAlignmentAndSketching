import sys
import os
from glob import glob
from os.path import join
import numpy as np

PROJ_DIR = '/home/gcgreen2/alignment'
GIT_DIR = join(PROJ_DIR,'SequenceAlignmentAndSketching')


sys.path.append(GIT_DIR)
import minhash, suffix_hash

DATASETS = 'test2  gc0.2_del3.16  gc0.2_sub1.68_ins8.04_del3.16  gc0.2  gc0.2_ins8.04  gc0.2_ins8.04_del3.16  gc0.2_sub1.68  gc0.2_sub25  malaria_reads  ecoli_k12'
DATASETS = DATASETS.split()
gt_path_ = lambda dset: join(GIT_DIR, 'data', dset, 'ground_truth.txt')
reads_path_ = lambda dset: join(GIT_DIR, 'data', dset, 'reads.fasta')
out_dir_ = lambda dset: join(PROJ_DIR, 'out', dset)
aln_path_ = lambda out_dir,k: join(out_dir, 'mh_k'+str(k)+'.tsv')

### PARAMS ##########################
date = '7-26'
METHOD = minhash
# args = {'methods':['inflec2','inflec','max','mld'], 'n_hash':500} 
args = {'ks':[i for i in range(10,20)], 'n_bits':24, 'n_hash':500} 
#####################################

for dset in DATASETS:   
    np.random.seed(0)
#     args['gt_path'] = gt_path_(dset)
    reads_path = reads_path_(dset)
    out_dir = out_dir_(dset)
    os.makedirs(out_dir, exist_ok=True)
    os.system(f'cp {join(GIT_DIR,"run_module_mh.py")} {join(out_dir,"run_module_mh_"+date+".py")}')
    aln_paths = [aln_path_(out_dir,k) for k in args['ks']]
    METHOD.find_overlaps(reads_path, aln_paths, **args)
    
   
