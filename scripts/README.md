# scripts

- Most scripts require Python >= 3.9 and pandas

## [make_votu_table.py](scripts/make_votu_table.py)

makes vOTU table from clustered vOTUs based on read mapping abundances.  Can output read counts, counts per million, or reads per kilobase

```bash
usage: make_votu_table.py [-h] [-v {count,cpm,rpk}] [-s SUFFIX] mmseqstsv filelist

positional arguments:
  mmseqstsv             mmseqs flattened tsv file (output of flatten_mmseqs_tsv.py)
  filelist              file containing list of sample abundance/count files (one per line) to include
                        in vOTU table. count files are output of verse. sample names should be
                        prefixes

options:
  -h, --help            show this help message and exit
  -v {count,cpm,rpk}, --value {count,cpm,rpk}
                        what value from abundance file to fill in the table (default: count)
  -s SUFFIX, --suffix SUFFIX
                        suffix of abundance files to remove when making table (default:
                        _virus.count.CDS.cpm.txt)
```



## [repseq_genomad_annotations.py](scripts/repseq_genomad_annotations.py)

makes table of taxonomic (and marker) gene annotations for representative sequences of vOTUs

```bash
usage: repseq_genomad_annotations.py [-h] votutable filelist

positional arguments:
  votutable   vOTU table
  filelist    file containing list of sample genomad virus_summary.tsv files (one per line) from which
              to get the repseq annotations. sample name should be prefix

options:
  -h, --help  show this help message and exit
```



## [genomad_virus2gff.py](scripts/genomad_virus2gff.py)

Make gff-type table (gtf) from genomad summary output

```
usage: genomad_virus2gff.py [-h] [-m HEADERMAP] virussummary

positional arguments:
  virussummary          genomad *_virus_summary.tsv file

options:
  -h, --help            show this help message and exit
  -m HEADERMAP, --headermap HEADERMAP
                        genomad *_assembly_headermap.txt file. optional used in case original contig
                        names were changed.
```



## [genomad_genes2gff.py](scripts/genomad_genes2gff.py)

Make gff-type table (gtf) from genomad *\*virus_genes.tsv* and optionally, *\*plasmid_genes.tsv* files.

```bash
  genesfile             genomad *_virus_genes.tsv file

options:
  -h, --help            show this help message and exit
  -m HEADERMAP, --headermap HEADERMAP
                        genomad *_assembly_headermap.txt file. optional used in case original contig
                        names were changed.
  -p PLASMIDGENESFILE, --plasmidgenesfile PLASMIDGENESFILE
                        genomad *_plasmid_genes.tsv file
```

## [plotnine_heatmap.py](scripts/plotnine_heatmap.py)

Uses the [plotnine](https://plotnine.readthedocs.io/en/stable/index.html) python library to make heatmap from csv table.

```bash
usage: plotnine_heatmap.py [-h] [-a ABUND] [-t TITLE] [-d DROPCOLS] inputfile outputfile

This script uses the python library plotnine to make a heatmap from a csv table. The y-axis column should be the first one in the table after any columns are dropped.

positional arguments:
  inputfile             file with table from which to make heatmap
  outputfile            filename of output pdf

optional arguments:
  -h, --help            show this help message and exit
  -a ABUND, --abund ABUND
                        abundance value in table. used for labeling heatmap. the default is generic. (default: abund)
  -t TITLE, --title TITLE
                        Title of heatmap; double quoted
  -d DROPCOLS, --dropcols DROPCOLS
                        Columns to drop/remove before plotting - separated by commas; double quoted
```

