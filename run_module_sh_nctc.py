import sys
import os
from glob import glob
from os.path import join
import numpy as np

PROJ_DIR = '/home/gcgreen2/alignment'
GIT_DIR = join(PROJ_DIR,'SequenceAlignmentAndSketching')
READS_DIR = join(PROJ_DIR, 'spectral_jaccard_similarity/filtered_fasta')
READS_PATHS = glob(READS_DIR+'/*.fasta')

sys.path.append(GIT_DIR)
import minhash, suffix_hash

# gt_path_ = lambda dset: join(GIT_DIR, 'data', dset, 'ground_truth.txt')
# reads_path_ = lambda dset: join(GIT_DIR, 'data', dset, 'reads.fasta')
out_dir_ = lambda dset: join(PROJ_DIR, 'out', dset)
aln_path_ = lambda out_dir,method: join(out_dir, method+'.tsv')

### PARAMS ##########################
date = '8-4'
METHOD = suffix_hash
# args = {'methods':['inflec','inflec2','max'], 'n_hash':500} 
args = {'methods':['max2'], 'n_hash':500, 'rc':False} 
# args = {'k':20, 'n_bits':24, 'n_hash':500} 
#####################################

for reads_path in READS_PATHS:   
    np.random.seed(0)
    dset = ''.join(reads_path.split('_')[:-1])
    dset = dset.split('/')[-1]
    out_dir = out_dir_(dset)
    aln_paths = [aln_path_(out_dir, method) for method in args['methods']]
    os.system(f'cp {join(GIT_DIR,"run_module_sh.py")} {join(out_dir,"run_module_sh_"+date+".py")}')
    args['gt_path'] = join(PROJ_DIR, 'spectral_jaccard_similarity', 'groundTruths', dset+'_daligner_ground_truth.txt')
#     assert not os.path.exists(aln_path)
    METHOD.find_overlaps(reads_path, aln_paths, **args)
    
   
