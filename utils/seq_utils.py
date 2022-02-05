import os
import pyfastx

BASES = {'A','T','C','G'}

def clean_seqs(seqs, replace_nonbase=''):
    for i in range(len(seqs)):
        seq = seqs[i]
        for b in BASES: seq = seq.replace(b.lower(), b)
        for char in set(seq).difference({'A','T','C','G'}):
            seq = seq.replace(char, replace_nonbase)
        seqs[i] = seq
    return seqs
    
def get_fasta(file):
    return pyfastx.Fasta(file)

def get_fastq(file):
    return pyfastx.Fastq(file)

def get_seqs(file, fa_fq='fa'):
    f = get_fasta(file) if fa_fq=='fa' else get_fastq(file)
    seqs = [s.seq for s in list(f)]
    seqs = clean_seqs(seqs)
    return seqs

def get_seq(file, seq_id=None, longest=False): # default is first seq
    f = get_fasta(file)
    return f.longest.seq if longest else f[seq_id].seq if seq_id else f[0].seq
    
def get_longest_seq_id(file): # seq id of first seq
    f = get_fasta(file)
    return f.longest.name

########### WRITING FILES

def write_fasta(file, seqs): # seqs should be dict
    with open(file, 'w') as fh:
        for seq_id in seqs:
            fh.write('>' + seq_id + "\n")
            fh.write(seqs[seq_id] + "\n")
