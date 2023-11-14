#!/usr/bin/env python3

"""
Our workflow runs bbmap first to deduplicate sequences.  This creates deduplicated fasta with sequence header names like dup1>dup2.
These sequences are then cluster by mmseqs2.  mmseqs createtsv makes a 2 column file - first column rep seq, second column seq (one row for each seq).
This script takes in as argument the output of mmseqs createtsv and returns a file in the same format with the deduped sequences each on their own line.

"""

import sys


infile = sys.argv[1]

rep_seq = ""

with open(infile, 'r') as f:
    for line in f:
        l = line.split()
        seqs = l[0].split(">")
        seqs = list(set(seqs + l[1].split(">")))
        rep = seqs[0]
        for s in seqs:
            print(rep, s, sep='\t')
