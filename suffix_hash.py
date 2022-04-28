import numpy as np
import random
import math
from scipy import stats
import string
import sys
from multiprocessing import Pool

sys.path.append('/home/gcgreen2/alignment/SequenceAlignmentAndSketching/utils')
import seq_utils as su
# import hash_utils as hu
import aln_utils as au


# The following script is used to evaluate the locational hashes of two DNA sequences.


#Function that converts string DNA alphabets (A,T,G,C) to numpy array numbers (0,1,2,3) respectively
def convert_DNA_to_numbers(string):
    temp = string.replace("A","0").replace("T","1").replace("G","2").replace("C","3")
    return np.array(list(temp), dtype=int)

# Function to Generate masks. The output is a 2D numpy array of masks
def obtain_masks(string_size,number_of_masks): #Select string_size to be the larger of the two strings
    masks = [np.random.choice([0,1],2*string_size) for i in range(number_of_masks)]
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

# def in_window(loc_diffs, diff, window_len): # diffs which are in the range [diff, diff+window+len]
#     return np.abs(loc_diffs-diff-window_len/2) <= window_len/2 

def get_mode(loc_diffs):#, window_len=100):
    quantiles = [round(diff,-1) for diff in loc_diffs] # round to nearest 10
    diffs, counts = np.unique(quantiles, return_counts=True); diffs,counts = list(diffs),list(counts)
    mode = diffs[np.argmax(counts)]
    c_low, c, c_high = [counts[diffs.index(mode+offset)] if mode+offset in diffs else 0 for offset in [-10,0,10]]
    mode_interpol = (c_low*(mode-10) + c*mode + c_high*(mode+10))//(c_low+c+c_high)
#     counts = [sum(in_window(loc_diffs, diff, window_len)) for diff in loc_diffs] # num diffs in each window
#     mode_diffs = loc_diffs[np.where(in_window(loc_diffs, loc_diffs[np.argmax(counts)], window_len))] # get the "clump" of diffs in the mode
    return mode_interpol, c_low+c+c_high
    
def est_overlap(sketches1, sketches2, len1, len2):
    loc_diffs = np.array([int(i-j) for i,j in zip(sketches1,sketches2)])
    mode, n_collisions = get_mode(loc_diffs)
    if (len1+len2)/(1+len(loc_diffs)/n_collisions) < 0.2*min(len1,len2): # number of collisions must imply at least 20% overlap
        return 0
    theta_hat = min(len1-mode, len2+mode) # one is the overlap, one is the fulll length of them combined
    return theta_hat 
    
# function to obtain sketches based on location hashing
def get_seq_sketches(seq, masks): #Sketch size is B
    seq = convert_DNA_to_numbers(seq);
    sketches = [lexicographic_first(seq,mask) for mask in masks]
    return sketches

def get_all_sketches(seqs, n_hash, n_bits):
    masks = obtain_masks(n_bits, n_hash);
    args = ((seq,masks) for seq in seqs)
    with Pool() as pool:
        sketches = pool.starmap(get_seq_sketches, args)
    return sketches

def find_overlaps(fasta_file, aln_file, n_hash=100, n_bits=None, k=None):
    np.random.seed(0)
    seqs = su.get_seqs(fasta_file); seq_lens = [len(s) for s in seqs]
    sketches = get_all_sketches(seqs, n_hash, max(seq_lens))
    au.find_overlaps(aln_file, sketches, seq_lens, est_overlap)

