#!/usr/bin/env python3

import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser()
parser.add_argument("mmseqstsv", help="mmseqs tsv file")
parser.add_argument("filelist", help="list of abundance files for each sample - output of verse. sample name should be prefix")
parser.add_argument("-v", "--value", help="what value to fill in the table", choices = ['count', 'cpm', 'rpk'], default = 'count')
parser.add_argument("-s", "--suffix", help="suffix of abundance files to remove when making table", default = '_virus.count.CDS.cpm.txt')
args = parser.parse_args()


def read_abund(abundfile, mmseqs, value, sfx=""):
    """abundfile: abundance filename - output of verse
       mmseqs: dataframe mapping repseq to seq
       value: which value column to include
       sfx: suffix to remove from abundfile for sample name
    """
    a = pd.read_csv(abundfile, sep='\t')
    a = a.merge(mmseqs, how='inner', left_on = 'gene', right_on = 'seq')
    a = a[['repseq', value]]
    samplename = abundfile.split("/")[-1]
    a['sample'] = samplename.removesuffix("_virus.count.CDS.cpm.txt")
    return(a)



mmseqs = pd.read_csv(args.mmseqstsv, sep='\t', names=['repseq', 'seq'], header=None)

abundlist = pd.read_csv(args.filelist, sep='\t', names = ['file'], header=None)
longtab = abundlist.groupby('file', group_keys=False).apply(lambda x: read_abund(abundfile = x.iloc[0]['file'], mmseqs = mmseqs, value = args.value, sfx = args.suffix))
votu = pd.pivot_table(longtab, values = args.value, index = 'repseq', columns = 'sample', fill_value = 0, aggfunc = "sum")
votu.to_csv(sys.stdout, sep = "\t")
