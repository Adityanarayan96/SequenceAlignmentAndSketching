import sys
import os
from glob import glob
from os.path import join

PROJ_DIR = '/home/gcgreen2/alignment'
GIT_DIR = join(PROJ_DIR,'SequenceAlignmentAndSketching')
DATASETS = 'gc0.2  gc0.2_del3.16  gc0.2_ins8.04  gc0.2_ins8.04_del3.16  gc0.2_sub1.68  gc0.2_sub1.68_ins8.04_del3.16  nctc1080  unif'
DATASETS = DATASETS.split()
#READS_DIR = join(PROJ_DIR, 'spectral_jaccard_similarity/filtered_fasta')
#READS_PATHS = glob(READS_DIR+'/*.fasta')

sys.path.append(GIT_DIR)
import minhash, loc_minhash

### PARAMS ##########################
k = 12
n_bits = 24
n_hash = 200
OUT_DIR = join(PROJ_DIR, 'out', '4-14-lmh')
METHOD = loc_minhash
#####################################

os.makedirs(OUT_DIR, exist_ok=True)
for dset in DATASETS:   
    reads_path = join(GIT_DIR, 'data', dset, 'reads.fasta')
    aln_path = join(OUT_DIR, dset+'_aln.tsv')
    METHOD.find_overlaps(reads_path, aln_path, k=k, n_hash=n_hash, n_bits=n_bits)
