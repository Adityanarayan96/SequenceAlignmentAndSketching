import numpy as np
from tqdm import tqdm
from multiprocessing import Pool
import time
import sys

sys.path.append('/home/gcgreen2/alignment/SequenceAlignmentAndSketching/utils')
import seq_utils as su
import hash_utils as hu
   
def minhash(seq, k, h): # returns (minhash, location, is_revcomp)
    hashes, hashes_rc = h(seq, k)
    return min(min(hashes), min(hashes_rc))

def get_all_minhashes(seqs, k, n_hash, n_bits):
    get_seq_minhashes = lambda seq,k,hs: [minhash(seq,k,h) for h in hs]
    hs = [hu.random_hash_func(n_bits=n_bits) for _ in range(n_hash)]
    args = ((seqs[i],k,hs) for i in range(len(seqs)))
    with Pool() as pool:
        minhashes = pool.starmap(get_seq_minhashes, args)
    return minhashes

def overlap_est(minhashes1, minhashes2, len1, len2):
    alpha = sum([i==j for i,j in zip(minhashes1,minhashes2)])
    alpha /= len(minhashes1)
    theta_hat = (len1+len2) * alpha/(1+alpha)
    theta_hat = round(theta_hat, 1) 
    return theta_hat 
   
def pairwise_overlap_ests(seqs, minhashes):
    pairwise_ests = []
    for i in tqdm(range(len(seqs)), desc='Estimating pairwise overlaps', leave=True):
        for j in range(len(seqs)):
            if i==j: continue
            theta_hat = overlap_est(minhashes[i], minhashes[j], len(seqs[i]), len(seqs[j]))
            if theta_hat > 0: 
                pairwise_ests.append((i+1,j+1,theta_hat,len(seqs[i]),len(seqs[j])))
    return pairwise_ests
    
def find_overlaps(fasta_file, aln_file, k, n_hash=100, n_bits=20):
    np.random.seed(0)
    seqs = su.get_seqs(fasta_file)
    minhashes = get_all_minhashes(seqs, k, n_hash, n_bits)
    pairwise_ests = pairwise_overlap_ests(seqs, minhashes)
    with open(aln_file, 'w') as fh:
        for aln in pairwise_ests:
            fh.write('\t'.join([str(x) for x in aln])+'\n')

