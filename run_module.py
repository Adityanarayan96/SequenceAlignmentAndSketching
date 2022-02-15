import sys
import os
from os.path import join

PROJ_DIR = '/home/gcgreen2/alignment'
OUT_DIR = join(PROJ_DIR, 'out', 'mh')
READS_DIR = join(PROJ_DIR, 'spectral_jaccard_similarity/filtered_fasta')

sys.path.append(join(PROJ_DIR,'SequenceAlignmentAndSketching'))
import minhash
os.makedirs(OUT_DIR, exist_ok=True)

# PARAMS
k = 16
n_bits = 24
n_hash = 200

for reads_file in os.listdir(READS_DIR):
    reads_path = join(READS_DIR, reads_file)
    aln_file = reads_file.replace('_reads.fasta','_aln.tsv')
    aln_path = join(OUT_DIR, aln_file)
    minhash.find_overlaps(reads_path, aln_path, k=k, n_bits=n_bits, n_hash=n_hash)