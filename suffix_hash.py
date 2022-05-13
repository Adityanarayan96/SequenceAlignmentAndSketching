import numpy as np
import random
import math
from scipy import stats
from tqdm import tqdm
import string
import sys
from multiprocessing import Pool

sys.path.append('/home/gcgreen2/alignment/SequenceAlignmentAndSketching/utils')
# import seq_utils as su
# import aln_utils as au


# The following script is used to evaluate the locational hashes of two DNA sequences.


#Function that converts string DNA alphabets (A,T,G,C) to numpy array numbers (0,1,2,3) respectively
def convert_DNA_to_numbers(string):
    temp = string.replace("A","0").replace("T","1").replace("G","2").replace("C","3")
    return np.array(list(temp), dtype=int)

# Function to Generate masks. The output is a 2D numpy array of masks
def obtain_masks(string_size,number_of_masks): #Select string_size to be the larger of the two strings
    masks = [np.random.choice([0,1,2,3],string_size) for i in range(number_of_masks)]
    return masks

def lexicographic_first(array,mask):
    best,idxs = [],[]
    for i,b in enumerate(array):
        idxs.append(i)
        best.append((array[i]-mask[i-idxs[0]]) % 4) # append to end of best (best idx is at idxs[0])
        j = 1
        while(len(idxs) > j):
            bm = (b-mask[i-idxs[j]]) % 4
            if bm > best[i-idxs[j]]:
                del idxs[j]
            elif bm < best[i-idxs[j]]:
                best = [*best[:i-idxs[j]], bm]
                idxs = idxs[j:]
                j = 1
            else: 
                j += 1
    return idxs[0]

######################################################

# def in_window(loc_diffs, diff, window_len): # diffs which are in the range [diff, diff+window+len]
#     return np.abs(loc_diffs-diff-window_len/2) <= window_len/2 

# def get_mode(loc_diffs):#, window_len=100):
#     quantiles = [round(diff,-1) for diff in loc_diffs] # round to nearest 10
#     diffs, counts = np.unique(quantiles, return_counts=True); diffs,counts = list(diffs),list(counts)
#     mode = diffs[np.argmax(counts)]
#     c_low, c, c_high = [counts[diffs.index(mode+offset)] if mode+offset in diffs else 0 for offset in [-10,0,10]]
#     mode_interpol = (c_low*(mode-10) + c*mode + c_high*(mode+10))//(c_low+c+c_high)
# #     counts = [sum(in_window(loc_diffs, diff, window_len)) for diff in loc_diffs] # num diffs in each window
# #     mode_diffs = loc_diffs[np.where(in_window(loc_diffs, loc_diffs[np.argmax(counts)], window_len))] # get the "clump" of diffs in the mode
#     return mode_interpol, c_low+c+c_high
    
# def est_overlap(sketches1, sketches2, len1, len2):
#     loc_diffs = np.array([int(i-j) for i,j in zip(sketches1,sketches2)])
#     mode, n_collisions = get_mode(loc_diffs)
#     if (len1+len2)/(1+len(loc_diffs)/n_collisions) < 0.2*min(len1,len2): # number of collisions must imply at least 20% overlap
#         return 0
#     theta_hat = min(len1-mode, len2+mode) # one is the overlap, one is the fulll length of them combined
#     return theta_hat 

###############################################

def edit_dist(s, t, d=5):
    if len(s)!=len(t): return -1
    rows = len(s)+1
#     deletes, inserts, substitutes = (1,1,1)
    
    dist = [[np.inf for _ in range(2*d+3)] for i in range(rows)]

    # init memoization
    for row in range(d+1):
        dist[row][d-row+1] = row * deletes
    for col in range(d+1):
        dist[0][col+d+1] = col * inserts 
        
    # dynamic programming
    for row in range(1, rows):
        for col in range(max(-d,-row+1),min(d+1,rows-row)):
            cost = 0 if s[row-1] == t[row+col-1] else 1
            dist[row][col+d+1] = min(dist[row-1][col+d+2] + 1,
                                 dist[row][col+d] + 1,
                                 dist[row-1][col+d+1] + cost) # substitution
#         if min(dist[row]) > 2*d: return -1
    return min(dist[-1])

def est_overlap_edit_dist(sketch1, sketch2, seq1, seq2, k=30):
    edit_dists = [edit_dist(seq1[i:i+k],seq2[j:j+k]) for i,j in zip(sketch1, sketch2)]
    similarity = sum([k-e for e in edit_dists if e!=-1])
    return similarity, edit_dists

def est_overlap_edit_dist_filter(sketch1, sketch2, seq1, seq2, k=30):
    edit_dists = [edit_dist(seq1[i:i+k],seq2[j:j+k]) for i,j in zip(sketch1, sketch2)]
    similarity = sum([k-e for e in edit_dists if e!=-1])
    return similarity, edit_dists

###################################################

def get_matching_bases(i, j, seq1, seq2):
    count = 0
    while(i<len(seq1) and j<len(seq2) and seq1[i]==seq2[j]):
        count += 1
        i += 1
        j += 1
    return count

def est_overlap_top_matching(sketch1, sketch2, seq1, seq2):
    n_matching = [get_matching_bases(i,j,seq1,seq2) for i,j in zip(sketch1, sketch2)]
    total_match = sum(n_matching)
    return total_match, n_matching

def est_overlap_top_matching_filter(sketch1, sketch2, seq1, seq2,c): #c is a parameter that increases the threshold as (1+c)*median
    n_matching = [get_matching_bases(i,j,seq1,seq2) for i,j in zip(sketch1, sketch2)]
    threshold = np.median(n_matching)*(1+c)
    n_matching = np.array(n_matching)
    n_matching[np.argwhere(n_matching <= threshold)] = 0
    n_matching = n_matching.tolist()
    total_match = sum(n_matching)
    return total_match, n_matching

# def est_overlap_top_matching_filter_with_sampling(sketch1, sketch2, seq1, seq2,c): #c is a parameter that increases the threshold as (1+c)*median
#     n_matching = [get_matching_bases(i,j,seq1,seq2) for i,j in zip(sketch1, sketch2)]
#     threshold = np.median(n_matching)*(1+c)
#     n_matching = np.array(n_matching)
#     n_matching[np.argwhere(n_matching <= threshold)] = 0
#     n_matching = n_matching.tolist()
#     total_match = sum(n_matching)
#     return total_match, n_matching 
######################################################
    
def pairwise_overlap_ests(sketches, seqs, seq_lens, est_overlap):
    estimates = []
    n = len(seqs)
    for i in tqdm(range(n-1), desc='Estimating pairwise overlaps', leave=True):
        for j in range(i+1,n):
            statistic, arr = est_overlap(sketches[i], sketches[j], seqs[i], seqs[j])
            if statistic > 0: 
                estimates.append((i+1,j+1,statistic,seq_lens[i],seq_lens[j],arr))
    return estimates

def write_overlaps(aln_file, estimates):
    print(f'# overlaps: {len(estimates)}')
    with open(aln_file, 'w') as fh:
        for aln in estimates:
            fh.write('\t'.join([str(x) for x in aln])+'\n')
            
######################################################
    
# function to obtain sketches based on location hashing
def get_seq_sketches(seq, masks): #Sketch size is B
    seq = convert_DNA_to_numbers(seq);
    sketches = [lexicographic_first(seq,mask) for mask in masks]
    return sketches

def get_all_sketches(seqs, n_hash): #Modified this to not include n_bits, since it is confusing
    seq_lens = [len(s) for s in seqs]
    masks = obtain_masks(max(seq_lens), n_hash);
    args = ((seq,masks) for seq in seqs)
    with Pool() as pool:
        sketches = pool.starmap(get_seq_sketches, args)
    return sketches

######################################################

def find_overlaps(fasta_file, aln_file, n_hash=100, n_bits=None, k=None, est_method=1):
    seqs = su.get_seqs(fasta_file); seq_lens = [len(s) for s in seqs]
    sketches = get_all_sketches(seqs, n_hash)
    if est_method==1: 
        threshold = get_threshold(sketches, seqs)
        estimates = pairwise_overlap_ests(sketches, seqs, seq_lens, est_overlap_top_matching)
    elif est_method==2:
        estimates = pairwise_overlap_ests(sketches, seqs, seq_lens, est_overlap_edit_dist)
        
    write_overlaps(aln_file, estimates)

