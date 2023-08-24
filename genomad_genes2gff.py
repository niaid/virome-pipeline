#!/usr/bin/env python3

"""
gff3 spec http://useast.ensembl.org/info/website/upload/gff3.html
gene vs CDS https://biology.stackexchange.com/questions/29929/what-is-the-difference-between-gene-and-cds-annotations

"""

import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("genesfile", help="genomad *_virus_genes.tsv file")
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


def make_gff(inputfile, plasmid=False):
    """
    Args:
      inputfile: genomad *_genes.tsv file
      plasmid: logical. is this a plasmid genes file? if so, Note=plasmid will be added to the gff attributes.

    Return:
      doesn't return anything.  prints to STDOUT
    """
    with open(inputfile, 'r') as gf:
        gfreader = csv.DictReader(gf, delimiter="\t")
        for row in gfreader:
            contig = row['gene'].rstrip('0123456789')
            contig = contig.removesuffix('_')
            # source = 'genomad'
            # type = 'CDS'
            # start = row['start']
            # end = row['end']
            score = '.' if row['bitscore'] == "NA" else row['bitscore']
            strand = strand_sign(row['strand'])
            # phase = '.'
            attributes = "ID=" + row['gene'] + ';Name=' + row['gene']
            if plasmid: attributes = attributes + ';Note=plasmid'
            print(contig,'genomad','CDS',row['start'],row['end'],score,strand,'.',attributes, sep='\t', flush=True)


## print header
print('##gff-version 3')

make_gff(args.genesfile)
if (args.plasmidgenesfile is not None):
    make_gff(args.plasmidgenesfile, plasmid=True)
