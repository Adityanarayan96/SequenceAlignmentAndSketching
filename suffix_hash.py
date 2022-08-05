import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm
from tqdm import tqdm
import pandas as pd
import string
from time import time
import sys
from multiprocessing import Pool

sys.path.append('/home/gcgreen2/alignment/SequenceAlignmentAndSketching/utils')
import seq_utils as su
import utils
# import aln_utils as au


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

def get_matching_bases(i, j, seq1, seq2):
    count = 0
    while(i<len(seq1) and j<len(seq2) and seq1[i]==seq2[j]):
        count += 1
        i += 1
        j += 1
    return count

def est_overlap_top_matching(sketch1, sketch2, seq1, seq2):
    n_matching = [get_matching_bases(i,j,seq1,seq2) for i,j in zip(sketch1, sketch2)]
#     total_match = sum([m for m in n_matching if m>=thresh])
    return n_matching

#######################################################

def get_thresh(sample_ks):
    all_ks = np.concatenate(sample_ks)
    unique_vals, counts = np.unique(all_ks, return_counts=True)
    thresh = [np.log(counts[i])-np.log(counts[i-1]) for i in range(1,len(counts))]
    thresh = [thresh[i]-thresh[i-1] for i in range(1,len(thresh))]
    thresh = [unique_vals[i+2] for i in range(1,len(thresh)) if np.sign(thresh[i])!=np.sign(thresh[i-1])][0]
    print('inflection point threshold: {}'.format(thresh))
    return thresh


def get_sample_ks(sketches, seqs, sample_size=10):
    sample_ks = []
    n = len(seqs)
    for i in np.random.choice(np.arange(n), size=sample_size, replace=False):
        for j in range(n):
            if i==j: continue
            arr = est_overlap_top_matching(sketches[i], sketches[j], seqs[i], seqs[j])
            sample_ks.append(arr)
    return sample_ks

def ks_to_stat(method, k_thresh, gt_path=None):
    if method == 'inflec':
        get_stat = lambda arr: sum([max(0,v-k_thresh) for v in arr])
    elif method == 'inflec2':
        get_stat = lambda arr: sum([v for v in arr if v>k_thresh])
    elif method == 'max':
        get_stat = lambda arr: max(arr)
    elif method == 'max2':
        get_stat = lambda arr: round(np.dot(np.sort(np.partition(arr,-5)[-5:])[::-1], [2**(-x) for x in range(5)]),1)
#     elif method == 'mld':
#         llr = get_llr(estimates, gt_path, ids)
#         get_stat = lambda arr: sum([llr[int(v)] for v in arr])
#     statistic = get_stat(ks)
    return get_stat

def get_stat_threshs(sample_ks, stat_funcs, quantile=0.75):
    method_stats = [[func(ks) for ks in sample_ks] for func in stat_funcs]
    stat_threshs = [np.quantile(stats, quantile) for stats in method_stats]
    return stat_threshs

def get_threshs_and_funcs(sketches, seqs, methods, rc):
    sample_ks = get_sample_ks(sketches, seqs) 
    k_thresh = get_thresh(sample_ks)
    stat_funcs = [ks_to_stat(method, k_thresh) for method in methods]
    stat_threshs = get_stat_threshs(sample_ks, stat_funcs, quantile=0.85 if rc else 0.7)
    return stat_threshs, stat_funcs

######################################################
    
def get_seq_sketches(seq, masks): 
    seq = convert_DNA_to_numbers(seq);
    sketches = [lexicographic_first(seq,mask) for mask in masks]
    return sketches

def get_all_sketches(seqs, masks):
    start = time()
    args = ((seq,masks) for seq in seqs)
    with Pool() as pool:
        sketches = pool.starmap(get_seq_sketches, args)
    print('sketching took {} minutes'.format(int((time()-start)/60)))
    return sketches

######################################################

def pairwise_overlap_ests(sketches, seqs, methods):
    n = len(seqs)
    method_ests = [[] for _ in methods]
    stat_threshs, stat_funcs = get_threshs_and_funcs(sketches, seqs, methods, rc=False)
    start = time()
    for i in range(n-1): 
        for j in range(i+1,n):
            ks = est_overlap_top_matching(sketches[i], sketches[j], seqs[i], seqs[j])
            for m in range(len(methods)):
                stat = stat_funcs[m](ks)
                if stat > stat_threshs[m]:
                    method_ests[m].append((i+1,j+1,stat,len(seqs[i]),len(seqs[j]),'+'))
                    method_ests[m].append((j+1,i+1,stat,len(seqs[j]),len(seqs[i]),'+'))
        if i==0: print('estimated total time for alignment: {} min.'.format(int(n*(time()-start)/60)))
    print('pairwise alignment took {} minutes'.format(int((time()-start)/60)))
    return method_ests

def pairwise_overlap_ests_rc(sketches, sketches_rc, seqs, seqs_rc, methods):
    n = len(seqs)
    method_ests = [[] for _ in methods]
    stat_threshs, stat_funcs = get_threshs_and_funcs(sketches, seqs, methods, rc=True)
    start = time()
    for i in range(n-1):
        for j in range(i+1,n):
            ks = est_overlap_top_matching(sketches[i], sketches[j], seqs[i], seqs[j])
            ks_rc = est_overlap_top_matching(sketches[i], sketches_rc[j], seqs[i], seqs_rc[j])
            for m,method in enumerate(methods):
                stat = stat_funcs[m](ks)
                if stat > stat_threshs[m]:
                    method_ests[m].append((i+1,j+1,stat,len(seqs[i]),len(seqs[j]),'+'))
                    method_ests[m].append((j+1,i+1,stat,len(seqs[j]),len(seqs[i]),'+'))
                stat_rc = stat_funcs[m](ks_rc)
                if stat_rc > stat_threshs[m]:
                    method_ests[m].append((i+1,j+1,stat_rc,len(seqs[i]),len(seqs[j]),'-'))
                    method_ests[m].append((j+1,i+1,stat_rc,len(seqs[j]),len(seqs[i]),'-'))
        if i==0: print('estimated total time for alignment: {} min.'.format(int(n*(time()-start)/60)))
    print('pairwise alignment took {} minutes'.format(int((time()-start)/60)))
    return method_ests

def find_overlaps(fasta_file, aln_paths, n_hash, rc=False, **args):
    seqs,ids = su.get_seqs(fasta_file, return_ids=True)
    masks = obtain_masks(max([len(seq) for seq in seqs]), n_hash);
    sketches = get_all_sketches(seqs, masks)
    if rc == True:
        seqs_rc = su.revcomp_seqs(seqs)
        sketches_rc = get_all_sketches(seqs_rc, masks)
        method_ests = pairwise_overlap_ests_rc(sketches, sketches_rc, seqs, seqs_rc, args['methods'])
    else:
        method_ests = pairwise_overlap_ests(sketches, seqs, args['methods'])
    for estimates,aln_path in zip(method_ests, aln_paths):
        utils.write_overlaps(aln_path, estimates, ids)
    
    
    
    
#####################################################
#####################################################
#####################################################

# EDIT_K = 20
# NORMAL1 = {'sgn':1, 'init':(7,2,2), 'bounds':[(4,20),(1,5),(0.5,0.9)]}
# NORMAL2 = {'sgn':-1, 'init':(13,2,2), 'bounds':[(6,19),(0.5,5),(0.5,0.9)]}

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

# def edit_dist(s, t, d=4, thresh=EDIT_K):
#     if len(s)!=len(t): return -1
    
#     n_matching = get_matching_bases(0,0,s,t)
#     if len(s)-n_matching<=d: return 0
#     s,t = s[n_matching:],t[n_matching:]
    
#     rows = len(s)+1 
#     dist = [[np.inf for _ in range(2*d+3)] for i in range(rows)]

#     # init memoization
#     for row in range(d+1):
#         dist[row][d-row+1] = row
#     for col in range(d+1):
#         dist[0][col+d+1] = col 
        
#     # dynamic programming
#     for row in range(1, rows):
#         for col in range(max(-d,-row+1),min(d+1,rows-row)):
#             cost = 0 if s[row-1] == t[row+col-1] else 1
#             dist[row][col+d+1] = min(dist[row-1][col+d+2] + 1,
#                                  dist[row][col+d] + 1,
#                                  dist[row-1][col+d+1] + cost) # substitution
#     return min(dist[-1])

# def est_overlap_edit_dist(sketch1, sketch2, seq1, seq2, k=EDIT_K, thresh=EDIT_K):
#     edit_dists = [edit_dist(seq1[i:i+k],seq2[j:j+k]) for i,j in zip(sketch1, sketch2)]
#     total_dist = sum([e for e in edit_dists if e!=-1 and e<=thresh])
#     return total_dist, edit_dists

###################################################


# def get_llr(estimates, gt_path, ids):
#     gt_df = pd.read_csv(gt_path, sep='\t', header=None, names=['i1','i2','overlap','l1','l2'])
#     gt_tuples = set([tuple(i) for i in gt_df[['i1','i2']].values])
#     err_vals = np.concatenate([e[2] for e in estimates if id_tuple(e,ids) not in gt_tuples])
#     correct_vals = np.concatenate([e[2] for e in estimates if id_tuple(e,ids) in gt_tuples])
#     max_val = max(max(err_vals),max(correct_vals))
    
#     vals,counts = np.unique(err_vals, return_counts=True)
#     err_pmf = {v:p/sum(counts) for v,p in zip(vals,counts)}
#     min_err_prob = min(err_pmf.values())
#     vals,counts = np.unique(correct_vals, return_counts=True)
#     corr_pmf = {v:p/sum(counts) for v,p in zip(vals,counts)}
#     llr = [np.log((corr_pmf[v] if v in corr_pmf else 0)/(err_pmf[v] if v in err_pmf else min_err_prob)) for v in range(max_val+1)]
    
#     # set all values after max equal to the max times linearly growing coefficient
#     maxidx = np.argmax(llr)
#     llr = [*llr[:maxidx], *[llr[maxidx]*i/maxidx for i in range(maxidx,len(llr))]]
# #     for i in range(maxidx, len(llr)):
# #         llr[i] = llr[maxidx]
#     return llr


# def get_stat_thresh(sample_ks, stat_funcs):
#     if method == 'inflec':
#         stat_thresh = 2
#     elif method == 'inflec2':
#         stat_thresh = int(2 * thresh + 2)
#     elif method == 'max':
#         stat_thresh = int(thresh+2)
#     elif method == 'max2':
#         stat_thresh = int(thresh*31/16+2)
#     return stat_thresh


# id_tuple = lambda e,ids: (ids[int(e[0])-1],ids[int(e[1])-1])
