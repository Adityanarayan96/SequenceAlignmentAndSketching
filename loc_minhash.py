import numpy as np
import sys

sys.path.append('/home/gcgreen2/alignment/SequenceAlignmentAndSketching/utils')
import seq_utils as su
import hash_utils as hu
import aln_utils as au

WINDOW_LEN = 200

def loc_minhash(seq, h): # returns (minhash, location)
    hashes = h(seq)
    loc = np.argmin(hashes)
    return hashes[loc], loc

# def get_loc_diffs(loc_minhashes1, loc_minhashes2, len1, len2):
#     loc_diffs,loc_diffs_rc = [],[]
#     loc_minhashes1 = loc_minhashes1[0]
#     loc_minhashes2,loc_minhashes2_rc = loc_minhashes2
#     for lmh1,lmh2 in zip(loc_minhashes1,loc_minhashes2):
#         if lmh1[0]==lmh2[0]:
#             loc_diffs.append(lmh1[1]-lmh2[1])
#     for lmh1,lmh2 in zip(loc_minhashes1,loc_minhashes2_rc):
#         if lmh1[0]==lmh2[0]:
#             loc_diffs_rc.append(lmh1[1]-lmh2[1])
#     loc_diffs = loc_diffs if len(loc_diffs)>len(loc_diffs_rc) else loc_diffs_rc
#     return np.array(loc_diffs)
def get_loc_diffs(loc_minhashes1, loc_minhashes2, len1, len2):
    loc_diffs = []
    for lmh1,lmh2 in zip(loc_minhashes1,loc_minhashes2):
        if lmh1[0]==lmh2[0]:
            loc_diffs.append(lmh1[1]-lmh2[1])
    return np.array(loc_diffs)

def weight_collisions(loc_diffs):
    counts = [sum(np.abs(loc_diffs-diff)<WINDOW_LEN) for diff in loc_diffs]
#     weight = sum([c for c in counts])
    return max(counts)

def est_overlap(loc_minhashes1, loc_minhashes2, len1, len2):
#     count_thresh = COUNT_THRESH * min(len1, len2)
    loc_diffs = get_loc_diffs(loc_minhashes1, loc_minhashes2, len1, len2)
    if len(loc_diffs) == 0: return 0
    weight = weight_collisions(loc_diffs)
    alpha = weight/len(loc_minhashes1)
    theta_hat = (len1+len2) * alpha/(1+alpha) 
    return round(theta_hat, 1) 

def find_overlaps(fasta_file, aln_file, k, n_hash=100, n_bits=20):
#     global COUNT_THRESH
#     COUNT_THRESH = 4**-k * n_hash
    seqs = su.get_seqs(fasta_file); seq_lens = [len(s) for s in seqs]
    loc_minhashes = hu.get_all_sketches(loc_minhash, seqs, k, n_hash, n_bits)
    au.find_overlaps(aln_file, loc_minhashes, seq_lens, est_overlap)
   


