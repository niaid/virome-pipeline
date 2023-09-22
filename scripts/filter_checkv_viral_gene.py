#!/usr/bin/env python3

import argparse
import sys
import os

lengththres=10000

parser = argparse.ArgumentParser(description="filtered fasta is printed to stdout")
parser.add_argument('indir', help='checkv directory that contains viruses.fna and contamination.tsv', type=str)
parser.add_argument('-l', '--length', help=f'minimum length for "manual check" default:{lengththres}', type=str)

args = parser.parse_args()

## input files
contam=os.path.join(args.indir, "contamination.tsv")
fasta=os.path.join(args.indir, "combined.fna")

## output file
# if args.outfile is None:
#     outfile=os.path.join(args.indir, outfile)
# else
#     outfile=args.outfile

## go through contam and look for contigs to keep
## https://www.protocols.io/view/viral-sequence-identification-sop-with-virsorter2-5qpvoyqebg4o/v3?step=4
keep = set()

with open(contam, 'r') as f:
    next(f)
    for line in f:
        l = line.split()
        ## Keep1
        if int(l[3]) > 0:
            keep.add(l[0])
            continue
        ## Keep2
        hg = int(l[4])
        if hg == 0:
            keep.add(l[0])
            continue
        if (hg == 1) and (int(l[1]) >= lengththres):
            keep.add(l[0])

## print filtered contigs
ctg2keep=0

with open(fasta, 'r') as f:
    for line in f:
        if line.startswith('>'):
            contig=line.strip()[1:]
            if contig in keep:
                sys.stdout.write(line)
                good=1
            else:
                good=0
            continue
        if good == 1:
            sys.stdout.write(line)
