#!/usr/bin/env python3

"""
Right now this script is slow.  Will speed up later.
"""

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("votutable", help="vOTU table")
parser.add_argument("filelist", help="file containing list of sample genomad virus_summary.tsv files (one per line) from which to get the repseq annotations. sample name should be prefix")
args = parser.parse_args()

def read_summary(infile, repseqs):
    """reads in summary file and prints repseq lines
    infile: *virus_summary.tsv file
    repseqs: set of rep seqs to look for

    """
    with open(infile, 'r') as f:
        for line in f:
            s = line.split()[0]
            if s in repseqs:
                print(line, end='')

## make set of repseqs
with open(args.votutable, 'r') as vt:
    samples = vt.readline().rstrip().split()
    samples = set(samples[1:])
    repseqs = set()
    for line in vt:
        rep = line.split()[0]
        repseqs.add(rep)

## get list of genomad virus summary files
with open(args.filelist, 'r') as fl:
    files = fl.read().splitlines()


## print header
with open(files[0], 'r') as f:
    print(f.readline(), flush=True, end="")

## go through summary files and print matching lines
[read_summary(x, repseqs) for x in files]
