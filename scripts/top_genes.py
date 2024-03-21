#!/usr/bin/env python3

"""
Requires python >= 3.9 and pandas
"""

import argparse
import pandas as pd
import sys

parser = argparse.ArgumentParser(description="This script takes in a gene abundance table and prints to STDOUT a subset of the table with top genes by mean abundance or prevalence and mean abundance across samples.",
    formatter_class=lambda prog: argparse.HelpFormatter(prog, width=90))

parser.add_argument("genetab", help="input gene abundance table")
parser.add_argument("-n", "--ntop", help="number of top genes to keep in output table. (default: %(default)s)", type=int, default=10)
parser.add_argument("-m", "--metric", help="how to choose top genes. by mean abundance across samples, by prevalence first and then mean, or by median. (default: %(default)s)", choices=['mean', 'prev', 'median'], default='mean')
parser.add_argument("-c", "--idcol", help="column name with gene id. default will use the first column")
parser.add_argument("-i", "--ignore", help="other columns besides idcol to ignore when choosing top genes. they will still be included in the output. column names should be separated by commas; double quoted")
parser.add_argument("-d", "--debug", help="print debug messages", action='store_true')
args = parser.parse_args()

## set debug cols
if args.debug: pd.set_option('display.max_columns', 500)

## get gene id column
if (args.idcol is None):
    args.idcol = 0

## read in gene table and set index to gene id column
## use python engine, so if there are blank columns they are skipped
df = pd.read_table(args.genetab, sep="\t+", index_col=args.idcol, engine='python')


if args.debug: print(f"dim of df {df.shape}", file=sys.stderr)


## drop rows that are all NaN
df.dropna(how='all', inplace=True)

## check number of rows if below number of top genes desired (or 0)
## we just return the same file
if (df.shape[0] < args.ntop):
    print(f"Number of genes in {args.genetab} is {df.shape[0]} which is below {args.ntop}. So printing entire table to STDOUT.", file=sys.stderr)
    df.to_csv(sys.stdout, sep = "\t")
    sys.exit(0)

## get cols to ignore and then columns to calculate mean/summary statistic
if (args.ignore is not None):
    toignore = [x.strip() for x in args.ignore.split(sep=',')]
    goodcols = list(set(df.columns) - set(toignore))
else:
    goodcols = df.columns

## calculate mean/median - and add as numbered column to avoid column name collision
rankcol = df.shape[1]
if (args.metric == 'median'):
    df[rankcol] = df[goodcols].median(axis=1, numeric_only=True).argsort()
else:
    df[rankcol] = df[goodcols].mean(axis=1, numeric_only=True).argsort()
sortby = [rankcol]

## calculate prevalence
if (args.metric == 'prev'):
    prevcol = rankcol + 1
    sortby = [prevcol] + sortby
    df[prevcol] = df[goodcols].gt(0).sum(axis=1)

if args.debug: print(f"columns to sort by {sortby}", file=sys.stderr)  

## get ntop rows and sort in descending order
df = df.sort_values(by=sortby, ascending=False).head(n=args.ntop)

if args.debug: 
    print("table with sort values", file=sys.stderr)
    print(df.head(), file=sys.stderr)

## remove mean/prev column
df = df.drop(labels = sortby, axis=1)

## write to STDOUT
df.to_csv(sys.stdout, sep = "\t")