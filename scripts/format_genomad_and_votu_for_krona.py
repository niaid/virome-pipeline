#!/usr/bin/env python3

import pandas as pd
import argparse
import os.path

parser = argparse.ArgumentParser(
                    description='Convert vOTU and genomad virus summary to Krona import format')
parser.add_argument('votu',type=str,
                    help='vOTU table')
parser.add_argument('summary',type=str,
                    help='path to repseq_genomad_virus_summary.tsv')
parser.add_argument('outdir',type=str,
                    help='path to output/temp dir for individual krona inputs')

args = parser.parse_args()

votu = pd.read_csv(args.votu, delimiter='\t', index_col = 'repseq')
tax = pd.read_csv(args.summary, delimiter='\t', index_col = 'seq_name', usecols=['seq_name', 'taxonomy'])

tax = tax['taxonomy'].str.split(';', expand=True)


def merge_and_write(v,t):
    s = str(v.columns[0])
    v = v.join(t)
    # subset to get nonzero rows
    v = v[v[s] > 0]
    v.to_csv(os.path.join(args.outdir, s + '.txt'), sep = '\t', index=False, header=False)

for sample in list(votu.columns):
    merge_and_write(votu[[sample]], tax)
