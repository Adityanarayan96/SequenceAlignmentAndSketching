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

### PARAMS ##########################
k = 20
n_bits = 26
n_hash = 500
OUT_DIR = join(PROJ_DIR, 'out', '6-20-sh-max')
METHOD = suffix_hash
args = {'method':'max'}  # estimation method 1: num front matching bases, 2: edit distance
#####################################

os.makedirs(OUT_DIR, exist_ok=True)
os.system(f'cp {join(GIT_DIR,"run_module2.py")} {OUT_DIR}')
for reads_path in READS_PATHS:   
    np.random.seed(0)
    aln_path = ''.join(reads_path.split('_')[:-1])
    aln_path = aln_path.split('/')[-1]
    aln_path = join(OUT_DIR, aln_path+'_aln.tsv')
    assert not os.path.exists(aln_path)
    METHOD.find_overlaps(reads_path, aln_path, k=k, n_hash=n_hash, n_bits=n_bits, **args)
    
   
