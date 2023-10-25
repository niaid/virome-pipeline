#!/usr/bin/env python3

"""
Requires python >= 3.9 and pandas

"""

import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser()
parser.add_argument("mmseqstsv", help="mmseqs flattened tsv file (output of flatten_mmseqs_tsv.py)")
parser.add_argument("filelist", help="file containing list of sample abundance/count files (one per line) to include in vOTU table. count files are output of verse. sample names should be prefixes")
parser.add_argument("-v", "--value", help="what value from abundance file to fill in the table (default: %(default)s)", choices = ['count', 'cpm', 'rpk'], default = 'count')
parser.add_argument("-s", "--suffix", help="suffix of abundance files to remove when making table (default: %(default)s)", default = '_virus.count.CDS.cpm.txt')
args = parser.parse_args()


def read_abund(abundfile, mmseqs, value, sfx=""):
    """abundfile: abundance filename - output of verse
       mmseqs: dataframe mapping repseq to seq
       value: which value column from abundfile to include
       sfx: suffix to remove from abundfile for sample name
    """
    a = pd.read_csv(abundfile, sep='\t')
    a = a.merge(mmseqs, how='inner', left_on = 'gene', right_on = 'seq')
    a = a[['repseq', value]]
    samplename = abundfile.split("/")[-1]
    a['sample'] = samplename.removesuffix(sfx)
    return(a)



mmseqs = pd.read_csv(args.mmseqstsv, sep='\t', names=['repseq', 'seq'], header=None)

abundlist = pd.read_csv(args.filelist, sep='\t', names = ['file'], header=None)
longtab = abundlist.groupby('file', group_keys=False).apply(lambda x: read_abund(abundfile = x.iloc[0]['file'], mmseqs = mmseqs, value = args.value, sfx = args.suffix))
votu = pd.pivot_table(longtab, values = args.value, index = 'repseq', columns = 'sample', fill_value = 0, aggfunc = "sum")
votu.to_csv(sys.stdout, sep = "\t")
