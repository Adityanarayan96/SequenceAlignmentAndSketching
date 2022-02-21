import numpy as np
from tqdm import tqdm
from multiprocessing import Pool
import time
import sys

sys.path.append('/home/gcgreen2/alignment/SequenceAlignmentAndSketching/utils')
import seq_utils as su
import hash_utils as hu

def loc_minhash(seq, k, h): # returns (minhash, location, is_revcomp)
    hashes, hashes_rc = h(seq, k)
    loc, loc_rc = [np.argmin(x) for x in [hashes,hashes_rc]]
    if hashes[loc] < hashes_rc[loc_rc]:
        return (hashes[loc], loc, False)
    else:
        return (hashes_rc[loc_rc], loc_rc, True)
    
def get_seq_minhashes(seq,k,hs): 
    return [loc_minhash(seq,k,h) for h in hs]

def get_all_loc_minhashes(seqs, k, n_hash, n_bits):
    hs = [hu.random_hash_func(n_bits=n_bits) for _ in range(n_hash)]
    args = ((seqs[i],k,hs) for i in range(len(seqs)))
    with Pool() as pool:
        loc_minhashes = pool.starmap(get_seq_minhashes, args)
    return loc_minhashes

def overlap_est(loc_minhashes1, loc_minhashes2, len1, len2):
    loc_diffs = [lmh1[1]-lmh2[1] for lmh1,lmh2 in zip(loc_minhashes1,loc_minhashes2) if lmh1[0]==lmh2[0]]
    return loc_diffs
#     alpha /= len(loc_minhashes1)
#     theta_hat = (len1+len2) * alpha/(1+alpha)
#     theta_hat = round(theta_hat, 1) 
#     return theta_hat 
   
def pairwise_overlap_ests(seqs, loc_minhashes):
    pairwise_ests = []
    for i in tqdm(range(len(seqs)), desc='Estimating pairwise overlaps', leave=True):
        for j in range(len(seqs)):
            if i==j: continue
            theta_hat = overlap_est(loc_minhashes[i], loc_minhashes[j], len(seqs[i]), len(seqs[j]))
            if theta_hat > 0: 
                pairwise_ests.append((i+1,j+1,theta_hat,len(seqs[i]),len(seqs[j])))
    return pairwise_ests
    
def find_overlaps(fasta_file, aln_file, k, n_hash=100, n_bits=20):
    np.random.seed(0)
    seqs = su.get_seqs(fasta_file)
    loc_minhashes = get_all_loc_minhashes(seqs, k, n_hash, n_bits)
    pairwise_ests = pairwise_overlap_ests(seqs, loc_minhashes,)
    with open(aln_file, 'w') as fh:
        for aln in pairwise_ests:
            fh.write('\t'.join([str(x) for x in aln])+'\n')

