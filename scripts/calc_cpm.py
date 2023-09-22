#!/usr/bin/env python3
"""
https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/https://www.rna-seqblog.com/rpkm-fpkm-and-tpm-clearly-explained/

takes in output of verse which has 3 columns - gene, count, length
and outputs to stdout 2 additional columns rpk and tpm

"""


import pandas as pd
import sys

df = pd.read_csv(sys.argv[1], sep='\t')
df['rpk'] = 1000 * df['count'] / df['length']
scalefactor = sum(df['rpk']) / 1e6
df['cpm'] = df['rpk'] / scalefactor

df.to_csv(sys.stdout, sep='\t', index=False)
