#!/usr/bin/env python3

"""
Requires python >= 3.9 and pandas and numpy
"""

import argparse
import pandas as pd
import numpy as np 
import sys

parser = argparse.ArgumentParser(
                    description='Create CPM table of hosts')
parser.add_argument('votu',type=str,
                    help='path to vOTU table')
parser.add_argument('iphop',type=str,
                    help='path to Host_prediction_to_genus_*.csv')

args = parser.parse_args()

##import vOTU table and host calls 
host = pd.read_csv(args.iphop, delimiter=',')
#host.to_csv("checking_host_file.csv")

votu = pd.read_csv(args.votu, delimiter='\t')

hostselect = host[['Virus', 'Host genus']]

## expand taxonomy calls 

hostselect[['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus']] = hostselect['Host genus'].str.split(';', expand=True)

#hostselect.to_csv("tax_split_check.tsv", sep="\t")

## If there are two genus calls for the same viral genome, replace genus calls with _g so the call will be made at the family level instead 

hostselect['Genus'] = np.where(hostselect['Virus'].duplicated(), 'g_', hostselect['Genus'])
#hostselect.to_csv("g_label_on_dups_check.tsv", sep="\t")

## merge data frames on genome name 
merged = votu.merge(hostselect, left_on='repseq', right_on='Virus')
#merged.to_csv("merged_votu_host_check.tsv", sep = "\t")

## remove duplicates (genomes with multiple genus calls) while keeping the sequence with the _g 
dedupe = merged.drop_duplicates('repseq', keep='last')

## Create updated  taxonomy 
dedupe['Host_genus'] = dedupe[['Domain','Phylum', 'Class', 'Order', 'Family', 'Genus']].agg(';'.join, axis=1)
#dedupe.to_csv("dropped_addcolumntest.tsv", sep = "\t")


# remove unnecessary columns that might be making an issue with sum
pruned = dedupe.drop([dedupe.columns[0], 'repseq', 'Virus', 'Host genus', 'Domain', 'Phylum',
'Class', 'Order', 'Family', 'Genus'], axis=1)


#pruned.to_csv("check_removal.tsv", sep = "\t")


# move the host_genus column to be the names of the rows 

pruned.insert(0, 'Host_genus', pruned.pop('Host_genus'))
#pruned.to_csv("pruned_reorder_test.tsv", sep = "\t")

## sum by host taxonomy to get CPM table of hosts 
sum = pruned.groupby(['Host_genus']).sum()
sum.to_csv(sys.stdout, sep = "\t")


