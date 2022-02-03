import os
import gzip
import pandas as pd
import numpy as np

BASES = {'A','T','C','G'}

def parse_fasta(genome):
    seqs = {}
    genome = genome.split('>')[1:]
    for seq in genome:
        seq = seq.split('\n')
        seq_id = seq[0].split(' ')[0]
        seq = ''.join(seq[1:])
        for char in set(seq).difference({'A','T','C','G'}):
            seq = seq.replace(char,'')
        seqs[seq_id] = seq
    return seqs
    
def get_seqs(file):
    if os.path.splitext(file)[1] == '.gz':
        with gzip.open(file) as fh:
            genome = fh.read().decode('utf-8').rstrip()
    else:
        with open(file) as fh:
            genome = fh.read().rstrip()
    seqs = parse_fasta(genome)
    return seqs

def get_seq(file, seq_id=None, largest=False): # default is first seq
    seqs = get_seqs(file)
    if seq_id is None:
        return list(seqs.values())[0]
    else:
        return seqs[seq_id]
    
def get_default_seq_id(genome_file): # seq id of first seq
    seqs = get_seqs(genome_file)
    return list(seqs.keys())[0]

def seq_len(file, seq_id=None):
    seq = get_seq(file, seq_id)
    return len(seq)

########### WRITING FILES

def write_fasta(file, seqs): # seqs should be dict
    with open(file, 'w') as fh:
        for seq_id in seqs:
            fh.write('>' + seq_id + "\n")
            fh.write(seqs[seq_id] + "\n")
