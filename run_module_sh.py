import sys
import os
from glob import glob
from os.path import join
import numpy as np

PROJ_DIR = '/home/gcgreen2/alignment'
GIT_DIR = join(PROJ_DIR,'SequenceAlignmentAndSketching')

sys.path.append(GIT_DIR)
import minhash, suffix_hash

# DATASETS = 'gc0.2_del3.16  gc0.2_sub1.68_ins8.04_del3.16  gc0.2  gc0.2_ins8.04  gc0.2_ins8.04_del3.16  gc0.2_sub1.68  gc0.2_sub25  malaria_sub1.68_ins8.04_del3.16'
DATASETS = 'malaria_reads  malaria_sub1.68_ins8.04_del3.16  ecoli_k12'
# DATASETS = 'test2'
DATASETS = DATASETS.split()
rc_dsets = ['malaria_reads','ecoli_k12','test2']
gt_path_ = lambda dset: join(GIT_DIR, 'data', dset, 'ground_truth.txt')
reads_path_ = lambda dset: join(GIT_DIR, 'data', dset, 'reads.fasta')
out_dir_ = lambda dset: join(PROJ_DIR, 'out', dset)
aln_path_ = lambda out_dir,method: join(out_dir, method+'.tsv')

### PARAMS ##########################
date = '8-2'
METHOD = suffix_hash
args = {'methods':['inflec','inflec2','max','max2'], 'n_hash':500} 
# args = {'methods':['max2'], 'n_hash':500} 
# args = {'k':20, 'n_bits':24, 'n_hash':500} 
#####################################

for dset in DATASETS: 
    print('running suffix hash on {}'.format(dset))
    np.random.seed(0)
    args['gt_path'] = gt_path_(dset)
    args['rc'] = dset in rc_dsets
    reads_path = reads_path_(dset)
    out_dir = out_dir_(dset)
    os.makedirs(out_dir, exist_ok=True)
    os.system(f'cp {join(GIT_DIR,"run_module_sh.py")} {join(out_dir,"run_module_sh_"+date+".py")}')
    aln_paths = [aln_path_(out_dir,method) for method in args['methods']]
    METHOD.find_overlaps(reads_path, aln_paths, **args)
    
   
