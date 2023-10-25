#!/usr/bin/env python3

"""
gff3 spec http://useast.ensembl.org/info/website/upload/gff3.html
gene vs CDS https://biology.stackexchange.com/questions/29929/what-is-the-difference-between-gene-and-cds-annotations

"""

import argparse
import csv
import re

parser = argparse.ArgumentParser()
parser.add_argument("genesfile", help="genomad *_virus_genes.tsv file")
parser.add_argument("-m", "--headermap", help="genomad *_assembly_headermap.txt file.  optional used in case original contig names were changed.")
parser.add_argument("-p", "--plasmidgenesfile", help="genomad *_plasmid_genes.tsv file", required=False)
args = parser.parse_args()


def strand_sign(n):
    """
    arg: strand value from genomad output
    return +/- for positive/negative strand as required by gff3
    """
    n = int(n)
    if (n<0):
        return('-')
    return('+')

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


def make_gff(inputfile, plasmid=False, headermap=None):
    """
    Args:
      inputfile: genomad *_genes.tsv file
      plasmid: logical. is this a plasmid genes file? if so, Note=plasmid will be added to the gff attributes.

    Return:
      doesn't return anything.  prints to STDOUT
    """
    if (headermap):
        hmap = read_headermap(headermap)
    with open(inputfile, 'r') as gf:
        gfreader = csv.DictReader(gf, delimiter="\t")
        for row in gfreader:
            contig = row['gene'].rstrip('0123456789')
            contig = contig.removesuffix('_')
            contig = re.sub("\\|provirus.+$", "", contig)
            if (headermap):
                contig = hmap[contig]
            # source = 'genomad'
            # type = 'CDS'
            # start = row['start']
            # end = row['end']
            ## start and end coordinates for provirus appear to be coordinates from the original contig
            score = '.' if row['bitscore'] == "NA" else row['bitscore']
            strand = strand_sign(row['strand'])
            # phase = '.'
            attributes = "ID=" + row['gene'] + ';Name=' + row['gene']
            if plasmid: attributes = attributes + ';Note=plasmid'
            print(contig,'genomad','CDS',row['start'],row['end'],score,strand,'.',attributes, sep='\t', flush=True)


## print header
print('##gff-version 3')

make_gff(args.genesfile, headermap=args.headermap)
if (args.plasmidgenesfile is not None):
    make_gff(args.plasmidgenesfile, plasmid=True, headermap=args.headermap)
