import numpy as np
from numpy.random import randint
from sympy import nextprime
from tqdm import tqdm
from multiprocessing import Pool
import time
import sys

sys.path.append('/home/gcgreen2/alignment/SequenceAlignmentAndSketching/utils')
import seq_utils as su

# N_THREADS = 50

# MinHash implementation

base_to_bin = {'A':'00', 'T':'01', 'G':'10', 'C':'11'}
rc_base = {'A':'T','T':'A','G':'C','C':'G'}
revcomp = lambda seq: ''.join([rc_base[b] for b in seq[::-1]])

class random_hash_func():
    def __init__(self,n_bits=20):
        self.a = randint(2**n_bits)
        self.b = randint(2**n_bits)
        self.max_hash = nextprime(2**n_bits)
        
    def __call__(self,bin_repr):
        hash_val = (int(bin_repr,2)*self.a + self.b) % self.max_hash
        return hash_val
            
def seq_to_bin(seq, rc=False):
    if rc: seq = revcomp(seq)
    bin_repr = ''.join([base_to_bin[b] for b in seq])
    return bin_repr

def minhash(seq, k, h):
    n_kmers = len(seq)-k+1; assert n_kmers>0
    bin_repr = seq_to_bin(seq)
    bin_repr_rc = seq_to_bin(seq, rc=True)
    hashes = [h(bin_repr[2*i:2*(i+k)]) for i in range(n_kmers)]
    hashes_rc = [h(bin_repr[2*i:2*(i+k)]) for i in range(n_kmers)]
    return min([*hashes, *hashes_rc])

def get_seq_minhashes(seq,k,hs):
    return [minhash(seq,k,h) for h in hs]

def get_all_minhashes(seqs, k, n_hash, n_bits):
    hs = [random_hash_func(n_bits=n_bits) for _ in range(n_hash)]
    args = ((seqs[i],k,hs) for i in range(len(seqs)))
    with Pool() as pool:
        minhashes = pool.starmap(get_seq_minhashes, args)
#     minhashes = [-1 for _ in range(len(seqs))]
#     n_rounds = int(np.ceil(len(seqs)/N_THREADS))
#     for r in tqdm(range(n_rounds), desc='Computing minhashes', leave=True):
#         threads = [threading.Thread(target=set_seq_minhashes, args=(minhashes,i,seqs[i],k,hs)) \
#                    for i in range(r*N_THREADS,min((r+1)*N_THREADS, len(seqs)))]
#         start=time.time()
#         for thread in threads: thread.start() 
#         print(time.time()-start)
#         for thread in threads: thread.join()
#         print(time.time()-start)
    return minhashes

def overlap_est(minhash1, minhash2, len1, len2, theta0):
    alpha = sum([i==j for i,j in zip(minhash1,minhash2)])
    alpha /= len(minhash1)
    theta_hat = (len1+len2) * alpha/(1+alpha)
    theta_hat = round(theta_hat, 1) #if theta_hat>theta0*min(len1,len2) else 0
    return theta_hat 
   
def pairwise_overlap_ests(seqs, minhashes, theta0):
    pairwise_ests = []
    for i in tqdm(range(len(seqs)), desc='Estimating pairwise overlaps', leave=True):
        for j in range(len(seqs)):
            if i==j: continue
            theta_hat = overlap_est(minhashes[i], minhashes[j], len(seqs[i]), len(seqs[j]), theta0)
            if theta_hat > 0: 
                pairwise_ests.append((i+1,j+1,theta_hat,len(seqs[i]),len(seqs[j])))
    return pairwise_ests
    
def find_overlaps(fasta_file, aln_file, k, n_hash=100, n_bits=20, theta0=0.1):
    np.random.seed(0)
    seqs = su.get_seqs(fasta_file)
    minhashes = get_all_minhashes(seqs, k, n_hash, n_bits)
    pairwise_ests = pairwise_overlap_ests(seqs, minhashes, theta0)
    with open(aln_file, 'w') as fh:
        for aln in pairwise_ests:
            fh.write('\t'.join([str(x) for x in aln])+'\n')

