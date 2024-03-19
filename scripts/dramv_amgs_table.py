#!/usr/bin/env python3

"""
Requires python >= 3.9 and pandas

"""

import argparse
import pandas as pd
import sys
import os.path

parser = argparse.ArgumentParser()
default_column = 'pfam_hits'
parser.add_argument("pipelinedir", help="directory containing pipeline outputs - including sample directories")
parser.add_argument("samplelist", help="file containing list of samples - should be names of sample directories that you want to include in the summary table - one per line")
parser.add_argument("-d", "--dramv", help="which dramv directory do you want to use the annotations from. e.g. dramv or amgs (default: %(default)s)", default="amgs")
parser.add_argument("-r", "--verse", help="which verse directory do you want to use the abundances from.  (default: %(default)s)", default="amgs/abund_amgs")
parser.add_argument("-v", "--value", help="what value from abundance file to fill in the table (default: %(default)s)", choices = ['count', 'cpm', 'rpk'], default = 'count')
parser.add_argument("-s", "--suffix", help="suffix for file with abundances (default: %(default)s)", default = "_amgs.count.gene.cpm.txt")
args = parser.parse_args()

print(args, file=sys.stderr)

def read_abund(abundfile, amgtable, value, sample):
    """abundfile: abundance filename - output of verse
       annotable: corresponding dramv annotations table from read_annotate
       value: which value column from abundfile to include
       sample: sample name associated with file
    """
    a = pd.read_csv(abundfile, sep='\t')
    a = a.merge(amgtable, how='inner', left_on = 'gene', right_on = 'gene')
    a = a[list(amgtable.columns) + [value]]
    return(a)


def read_amg(amgfile, sample):
    a = pd.read_csv(amgfile, sep='\t')
    ## get just the relevant columns
    column = ['gene', 'gene_id', 'gene_description']
    a = a[column]
    ## add sample column
    a['sample'] = sample
    return(a)

## read in samples
with open(args.samplelist, 'r') as f:
    samples = f.read().splitlines()

amgdf = read_amg(os.path.join(args.pipelinedir, samples[0], args.dramv, "dramv-distill", "amg_summary.tsv"), samples[0])


## make list of sample abundance files and amg summary files
abundlist = {sample: os.path.join(args.pipelinedir, sample, args.verse, sample + args.suffix) for sample in samples}
amglist = {sample: os.path.join(args.pipelinedir, sample, args.dramv, "dramv-distill", "amg_summary.tsv") for sample in samples}

## read in amg and abund files and make individual sample tables
amgtabs = {sample: read_amg(amglist[sample], sample) for sample in samples}
abundtabs = [read_abund(abundlist[sample], amgtabs[sample], args.value, sample) for sample in samples]

## concatenate abundance tables
longtab = pd.concat(abundtabs)

## pivot/cast table into wide gene_id table
finaltab = pd.pivot_table(longtab, values = args.value, index = ['gene_id', 'gene_description'], columns = 'sample', fill_value=0, aggfunc='sum')
finaltab.to_csv(sys.stdout, sep = "\t")
