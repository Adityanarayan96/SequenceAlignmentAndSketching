import numpy as np
import sys

sys.path.append('/home/gcgreen2/alignment/SequenceAlignmentAndSketching/utils')
import seq_utils as su
import hash_utils as hu
import aln_utils as au
   
def minhash(seq, h): # returns (minhash, location, is_revcomp)
    return h(seq)

def est_overlap(minhashes1, minhashes2, len1, len2):
    minhashes1 = minhashes1[0]
    minhashes2,minhashes2_rc = minhashes2
    alpha = sum([i==j for i,j in zip(minhashes1,minhashes2)])
    alpha_rc = sum([i==j for i,j in zip(minhashes1,minhashes2_rc)])
    alpha = max(alpha, alpha_rc)
    if alpha==0: return 0
    alpha /= len(minhashes1)
    theta_hat = (len1+len2) * alpha/(1+alpha)
    theta_hat = round(theta_hat, 1) 
    return theta_hat 

def find_overlaps(fasta_file, aln_file, k, n_hash=100, n_bits=20):
    seqs = su.get_seqs(fasta_file); seq_lens = [len(s) for s in seqs]
    minhashes = hu.get_all_sketches(minhash, seqs, k, n_hash, n_bits)
    au.find_overlaps(aln_file, minhashes, seq_lens, est_overlap)

