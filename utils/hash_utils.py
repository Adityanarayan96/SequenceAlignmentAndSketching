import os
import numpy as np
from numpy.random import randint
from sympy import nextprime

base_to_bin = {'A':'00', 'T':'01', 'G':'10', 'C':'11'}
rc_base = {'A':'T','T':'A','G':'C','C':'G'}
revcomp = lambda seq: ''.join([rc_base[b] for b in seq[::-1]])

class random_hash_func():
    def __init__(self,n_bits=20):
        self.a = randint(2**n_bits)
        self.b = randint(2**n_bits)
        self.max_hash = nextprime(2**n_bits)
    
    def seq_to_bin(self, seq, rc=False):
        if rc: seq = revcomp(seq)
        bin_repr = ''.join([base_to_bin[b] for b in seq])
        return bin_repr 
    
    def hash_val(self, bin_repr):
        hash_val = (int(bin_repr,2)*self.a + self.b) % self.max_hash
        return hash_val
    
    def __call__(self, seq, k):
        n_kmers = len(seq)-k+1; assert n_kmers>0
        bin_repr = self.seq_to_bin(seq)
        bin_repr_rc = self.seq_to_bin(seq, rc=True)
        hashes = [self.hash_val(bin_repr[2*i:2*(i+k)]) for i in range(n_kmers)]
        hashes_rc = [self.hash_val(bin_repr_rc[2*i:2*(i+k)]) for i in range(n_kmers)]
        return hashes, hashes_rc

