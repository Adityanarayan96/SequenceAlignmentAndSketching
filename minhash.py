import numpy as np
from sympy import nextprime
import sys
from multiprocessing import Pool

sys.path.append('/home/gcgreen2/alignment/SequenceAlignmentAndSketching/utils')
import seq_utils as su
import utils

base_to_bin = {'A':'00', 'T':'01', 'G':'10', 'C':'11'}
# rc_base = {'A':'T','T':'A','G':'C','C':'G'}
# revcomp = lambda seq: ''.join([rc_base[b] for b in seq[::-1]])

class random_hash_func():
    def __init__(self,n_bits,k):
        self.a = np.random.randint(2**n_bits)
        self.b = np.random.randint(2**n_bits)
        self.max_hash = nextprime(2**n_bits)
        self.k = k
    
    def seq_to_bin(self, seq):
#         if rc: seq = revcomp(seq)
        bin_repr = ''.join([base_to_bin[b] for b in seq])
        return bin_repr 
    
    def hash_val(self, bin_repr):
        hash_val = (int(bin_repr,2)*self.a + self.b) % self.max_hash
        return hash_val
    
    def __call__(self, seq):
        n_kmers = len(seq)-self.k+1; assert n_kmers>0
        bin_repr = self.seq_to_bin(seq)
        hashes = [self.hash_val(bin_repr[2*i:2*(i+self.k)]) for i in range(n_kmers)]
        return hashes
    
def get_seq_sketches(sketcher,seq,hs): 
    return [sketcher(seq,h) for h in hs]
#     return [sketcher(seq,h) for h in hs], [sketcher(revcomp(seq),h) for h in hs]

def get_all_sketches(seqs, k, n_hash, n_bits):
    hs = [random_hash_func(n_bits=n_bits, k=k) for _ in range(n_hash)]
    args = ((minhash,seqs[i],hs) for i in range(len(seqs)))
    with Pool() as pool:
        sketches = pool.starmap(get_seq_sketches, args)
    return sketches    
    
def minhash(seq, h): # returns (minhash, location, is_revcomp)
    return min(h(seq))

def est_overlap(minhashes1, minhashes2, seq1, seq2):
    alpha = sum([i==j for i,j in zip(minhashes1,minhashes2)])
    if alpha==0: return 0
    alpha /= len(minhashes1)
    theta_hat = (len(seq1)+len(seq2)) * alpha/(1+alpha)
    theta_hat = round(theta_hat, 1) 
    return theta_hat

def pairwise_overlap_ests(sketches, seqs):
    estimates = []
    n = len(seqs)
    for i in tqdm(range(n-1), desc='Estimating pairwise overlaps', leave=True):
        for j in range(i+1,n):
            est = est_overlap(sketches[i], sketches[j], seqs[i], seqs[j]) 
            estimates.append((i+1,j+1,est,len(seqs[i]),len(seqs[j]),'+'))
            estimates.append((j+1,i+1,est,len(seqs[j]),len(seqs[i]),'+'))
    return estimates

def pairwise_overlap_ests_rc(sketches, sketches_rc, seqs, seqs_rc):
    estimates = []
    n = len(seqs)
    for i in tqdm(range(n-1), desc='Estimating pairwise overlaps', leave=True):
        for j in range(i+1,n):
            est = est_overlap(sketches[i], sketches[j], seqs[i], seqs[j]) 
            estimates.append((i+1,j+1,est,len(seqs[i]),len(seqs[j]),'+'))
            estimates.append((j+1,i+1,est,len(seqs[j]),len(seqs[i]),'+'))
            
            est_rc = est_overlap(sketches[i], sketches_rc[j], seqs[i], seqs_rc[j]) 
            estimates.append((i+1,j+1,est_rc,len(seqs[i]),len(seqs[j]),'-'))
            estimates.append((j+1,i+1,est_rc,len(seqs[j]),len(seqs[i]),'-'))
    return estimates


def find_overlaps(fasta_file, aln_paths, ks, n_hash, n_bits, rc=False, **args):
    seqs,ids = su.get_seqs(fasta_file, return_ids=True)
    for aln_path,k in zip(aln_paths, ks):
        sketches = get_all_sketches(seqs, k, n_hash, n_bits)
        if rc == True:
            seqs_rc = su.revcomp_seqs(seqs)
            sketches_rc = get_all_sketches(seqs_rc, k, n_hash, n_bits)
            estimates = pairwise_overlap_ests_rc(sketches, sketches_rc, seqs, seqs_rc)
        else:
            estimates = pairwise_overlap_ests(sketches, seqs)
        utils.write_overlaps(aln_path, estimates, ids)
    
    
#######################################################
#ARCHIVE
    
# def pairwise_overlap_ests(seq_lens, sketches):
#     pairwise_ests = []
#     n = len(seq_lens)
#     for i in tqdm(range(n), desc='Estimating pairwise overlaps', leave=True):
#         for j in range(n):
#             if i==j: continue
#             theta_hat = est_overlap(sketches[i], sketches[j], seq_lens[i], seq_lens[j])
#             if theta_hat > 0: 
#     pairwise_ests.append((i+1,j+1,theta_hat,seq_lens[i],seq_lens[j]))
#     return pairwise_ests

# def find_overlaps(fasta_file, aln_file, k, n_hash=100, n_bits=20, **args):
#     seqs = su.get_seqs(fasta_file); seq_lens = [len(s) for s in seqs]
#     minhashes = hu.get_all_sketches(minhash, seqs, k, n_hash, n_bits)
#     au.find_overlaps(aln_file, minhashes, seq_lens, est_overlap)

