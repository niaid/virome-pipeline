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
parser.add_argument("-d", "--dramv", help="which dramv directory do you want to use the annotations from. e.g. dramv or amgs (default: %(default)s)", default="dramv")
parser.add_argument("-c", "--column", help=f"columns in dramv annotations file from which to make the gene table. can use more than once for more than one column (default: {default_column})", nargs='*', action='append')
parser.add_argument("-v", "--value", help="what value from abundance file to fill in the table (default: %(default)s)", choices = ['count', 'cpm', 'rpk'], default = 'count')
args = parser.parse_args()

## deal with column arg
if args.column is None:
    args.column = [default_column]
else:
    args.column = [item for row in args.column for item in row]

print(args, file=sys.stderr)

def read_abund(abundfile, annotable, value, sample):
    """abundfile: abundance filename - output of verse
       annotable: corresponding dramv annotations table from read_annotate
       value: which value column from abundfile to include
       sample: sample name associated with file
    """
    a = pd.read_csv(abundfile, sep='\t')
    a = a.merge(annotable, how='inner', left_on = 'gene', right_on = 'gene')
    a = a[list(annotable.columns) + [value]]
    return(a)

def read_annotate(annofile, column, sample):
    a = pd.read_csv(annofile, sep='\t')
    a = a.rename(columns={"Unnamed: 0": "gene"})
    ## remove rows where column values are empty https://stackoverflow.com/a/46764265
    a = a[(a[column].values == '').sum(axis=1) < len(column)]
    a = a.dropna(subset=column, thresh=len(column))
    ## get just the relevant columns
    column = ['gene'] + column
    a = a[column]
    ## add sample column
    a['sample'] = sample
    return(a)

## read in samples
with open(args.samplelist, 'r') as f:
    samples = f.read().splitlines()

## make list of sample files and dramv annotation files
abundlist = {sample: os.path.join(args.pipelinedir, sample, "verse_dramv", sample + "_genes.count.CDS.cpm.txt") for sample in samples}
annolist = {sample: os.path.join(args.pipelinedir, sample, args.dramv, "dramv-annotate", "annotations.tsv") for sample in samples}

annotabs = {sample: read_annotate(annolist[sample], args.column, sample) for sample in samples}

abundtabs = [read_abund(abundlist[sample], annotabs[sample], args.value, sample) for sample in samples]

longtab = pd.concat(abundtabs)

# print(longtab.tail(), file=sys.stderr)

finaltab = pd.pivot_table(longtab, values = args.value, index = args.column, columns = 'sample', fill_value=0, aggfunc='sum')
finaltab.to_csv(sys.stdout, sep = "\t")
