#!/usr/bin/env python3

"""
Uses pandas and plotnine v 0.12.4 python libraries
https://plotnine.readthedocs.io/en/stable/index.html


"""

import argparse
import pandas as pd
import sys

from plotnine import (
    ggplot,
    aes,
    geom_tile,
    geom_text,
    scale_y_reverse,
    scale_y_discrete,
    scale_fill_brewer,
    scale_color_manual,
    coord_equal,
    theme,
    theme_void,
    element_blank,
    element_rect,
    element_text,
    scale_fill_continuous,
    ggtitle
)
from mizani.transforms import log1p_trans


## arguments
parser = argparse.ArgumentParser(description="This script uses the python library plotnine to make a heatmap from a csv table.  The y-axis column should be the first one in the table after any columns are dropped.")
parser.add_argument("inputfile", help="file with table from which to make heatmap")
parser.add_argument("outputprefix", help="file prefix of output file(s)")
parser.add_argument("-a", "--abund", help="abundance value in table. used for labeling heatmap.  the default is generic. (default: %(default)s)", default="abund")
parser.add_argument("-t", "--title", help="Title of heatmap; double quoted")
parser.add_argument("-d", "--dropcols", help="Columns to drop/remove before plotting - separated by commas; double quoted")
parser.add_argument("-f", "--plotfiletype", help="File type of output plot. (default: %(default)s)", choices=['pdf', 'png', 'both'], default='both')
args = parser.parse_args()

print(args, file=sys.stderr)

## Load the CSV file
inputdf = pd.read_csv(args.inputfile, sep='\t', index_col=False)

## check if inputdf is empty
if inputdf.size == 0:
    print(f"WARNING: abundance file {args.inputfile} is empty. Heatmap will not be made.", file=sys.stderr)
    sys.exit(0)

if (args.dropcols is not None):
    todrop = args.dropcols.split(sep=',')
    inputdf = inputdf.drop(todrop, axis=1)

## first column is y-axis
y_axis = inputdf.columns[0]
inputdf = inputdf.set_index(y_axis)

## dimensions of plot
width = float(max(inputdf.shape[1] / 4, 5))
height = float(max((inputdf.shape[0] / 6.25) + 0.5, 5))


## melt dataframe
df = inputdf.stack().reset_index(name=args.abund)
df = df.rename(columns={'level_1': 'sample'})

## make plot
ggp = ggplot(df, aes('sample', y_axis, fill=args.abund)) + geom_tile(aes(width=.95, height=.95)) + coord_equal() + scale_fill_continuous(trans = log1p_trans) + theme(axis_text_x=element_text(rotation=90), figure_size=(width, height))  

## add title if needed
if (args.title is not None):
    ggp = ggp + ggtitle(title = args.title)


## save to file
if (args.plotfiletype in ['pdf', 'both']):
    ggp.save(filename=args.outputprefix + '.pdf', limitsize=False)

if (args.plotfiletype in ['png', 'both']):
    ggp.save(filename=args.outputprefix + '.png', limitsize=False)
