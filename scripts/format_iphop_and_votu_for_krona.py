#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np 
import os.path


parser = argparse.ArgumentParser(
                    description='Convert vOTU and iphop host genus prediction to Krona import format')
parser.add_argument('votu',type=str,
                    help='path to vOTU table')
parser.add_argument('iphop',type=str,
                    help='path to Host_prediction_to_genus_*.csv')
parser.add_argument('outdir',type=str,
                    help='path to output/temp dir for individual krona inputs')

args = parser.parse_args()

votu = pd.read_csv(args.votu, delimiter='\t', index_col = 'repseq')
tax = pd.read_csv(args.iphop, delimiter=',', usecols=['Virus', 'Host genus'])

#tax = tax['Host genus'].str.split(';', expand=True)
#tax.to_csv('outfile.csv', index=True)

tax[['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus']] = tax['Host genus'].str.split(';', expand=True)
#tax.to_csv('check_for_split.csv', index=True)

tax['Genus'] = np.where(tax['Virus'].duplicated(), 'g_', tax['Genus'])
tax.to_csv('check_for_replace.csv', index=True)

tax = tax.drop_duplicates('Virus', keep='last')
tax.to_csv('check_for_dedupe.csv', index=True)

tax.set_index('Virus', inplace=True)

def merge_and_write(v,t):
    s = str(v.columns[0])
    # keep only seqs in iphop taxonomy with how = 'right'
    v = v.join(t, how = 'right')
    # subset to get nonzero rows
    v = v[v[s] > 0]
    v.to_csv(os.path.join(args.outdir, s + '.txt'), sep = '\t', index=False, header=False)

for sample in list(votu.columns):
    merge_and_write(votu[[sample]], tax)
