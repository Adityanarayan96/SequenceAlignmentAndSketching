import sys
import os
from glob import glob
from os.path import join

PROJ_DIR = '/home/gcgreen2/alignment'
GIT_DIR = join(PROJ_DIR,'SequenceAlignmentAndSketching')
DATASETS = 'test'#'gc0.2  gc0.2_del3.16  gc0.2_ins8.04  gc0.2_ins8.04_del3.16  gc0.2_sub1.68  gc0.2_sub1.68_ins8.04_del3.16  nctc1080  unif  gc0.2_sub25  gc0.2_sub50'
DATASETS = DATASETS.split()
#READS_DIR = join(PROJ_DIR, 'spectral_jaccard_similarity/filtered_fasta')
#READS_PATHS = glob(READS_DIR+'/*.fasta')

sys.path.append(GIT_DIR)
import minhash, loc_minhash, suffix_hash

### PARAMS ##########################
k = 16
n_bits = 24
n_hash = 100
OUT_DIR = join(PROJ_DIR, 'out', '4-27-sh')
METHOD = suffix_hash
#####################################

os.makedirs(OUT_DIR, exist_ok=True)
os.system(f'cp {join(GIT_DIR,"run_module.py")} {OUT_DIR}')
for dset in DATASETS:   
    reads_path = join(GIT_DIR, 'data', dset, 'reads.fasta')
    aln_path = join(OUT_DIR, dset+'_aln.tsv')
    assert not os.path.exists(aln_path)
    METHOD.find_overlaps(reads_path, aln_path, k=k, n_hash=n_hash, n_bits=n_bits)
