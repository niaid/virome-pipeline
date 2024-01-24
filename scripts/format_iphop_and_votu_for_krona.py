#!/usr/bin/env python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser(
                    description='Convert vOTU and iphop host genus prediction to Krona import format')
parser.add_argument('votu',type=str,
                    help='path to vOTU table')
parser.add_argument('iphop',type=str,
                    help='path to Host_prediction_to_genus_*.csv')

args = parser.parse_args()

votu = pd.read_csv(args.votu, delimiter='\t', index_col = 'repseq')
tax = pd.read_csv(args.iphop, delimiter=',', index_col = 'Virus', usecols=['Virus', 'Host genus'])

tax = tax['Host genus'].str.split(';', expand=True)

# tax.to_csv('outfile.csv', index=True)

def merge_and_write(v,t):
    s = str(v.columns[0])
    # keep only seqs in iphop taxonomy with how = 'right'
    v = v.join(t, how = 'right')
    # subset to get nonzero rows
    v = v[v[s] > 0]
    v.to_csv(s + '.to_krona.txt', sep = '\t', index=False, header=False)

for sample in list(votu.columns):
    merge_and_write(votu[[sample]], tax)
