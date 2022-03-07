import sys
import os
from glob import glob
from os.path import join

PROJ_DIR = '/home/gcgreen2/alignment'
READS_DIR = join(PROJ_DIR, 'spectral_jaccard_similarity/filtered_fasta')
READS_PATHS = glob(READS_DIR+'/*.fasta')

sys.path.append(join(PROJ_DIR,'SequenceAlignmentAndSketching'))
import minhash
import loc_minhash

### PARAMS ##########################
k = 12
n_bits = 24
n_hash = 200
OUT_DIR = join(PROJ_DIR, 'out', f'mh_k{str(k)}_h{str(n_hash)}')
METHOD = minhash
#####################################

os.makedirs(OUT_DIR, exist_ok=True)
for reads_path in READS_PATHS:
    reads_file = os.path.basename(reads_path)
    aln_file = reads_file.replace('_reads.fasta','_aln.tsv')
    aln_path = join(OUT_DIR, aln_file)
    METHOD.find_overlaps(reads_path, aln_path, k=k, n_hash=n_hash, n_bits=n_bits)
