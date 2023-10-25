#!/usr/bin/env python3

"""
gff3 spec http://useast.ensembl.org/info/website/upload/gff3.html
gene vs CDS https://biology.stackexchange.com/questions/29929/what-is-the-difference-between-gene-and-cds-annotations

"""

import argparse
import csv
import re

parser = argparse.ArgumentParser()
parser.add_argument("virussummary", help="genomad *_virus_summary.tsv file")
parser.add_argument("-m", "--headermap", help="genomad *_assembly_headermap.txt file.  optional used in case original contig names were changed.")
args = parser.parse_args()



def read_headermap(headermap):
    hmap = {}
    with open(headermap, 'r') as f:
        for line in f:
            line = line.rstrip()
            l = line.split("\t")
            old = l[0].split()[0]
            new = l[1].split()[0]
            hmap[new] = old
    return(hmap)


def make_gff(inputfile, headermap=None):
    """
    Args:
      inputfile: genomad *_virus_summary.tsv file
      headermap: filename. file mapping from old contig names to new ones. optional

    Return:
      doesn't return anything.  prints to STDOUT
    """
    if (headermap):
        hmap = read_headermap(headermap)
    with open(inputfile, 'r') as gf:
        gfreader = csv.DictReader(gf, delimiter="\t")
        for row in gfreader:
            contig = re.sub("\\|provirus.+$", "", row['seq_name'])
            if (headermap):
                contig = hmap[contig]
            # source = 'genomad'
            # type = 'CDS'
            if row['coordinates'] == "NA":
                start = "1"
                end = row['length']
            else:
                [start, end] = row['coordinates'].split("-")
            score = row['virus_score']
            # strand = '+'
            # phase = '.'
            attributes = "ID=" + row['seq_name'] + ';Name=' + row['seq_name']
            print(contig,'genomad','CDS',start,end,score,'+','.',attributes, sep='\t', flush=True)


## print header
print('##gff-version 3')

make_gff(args.virussummary, headermap=args.headermap)
