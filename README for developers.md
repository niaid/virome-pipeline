# DiscoVir

## Files for running pipeline on LOCUS

-   Snakefile: pipeline script (reads in configs, commands for each
    pipeline step/rule)

-   cluster_setup.smk: helper script for reading in cluster config file

-   project_config.yaml: for snakemake --configfile option. config file
    with details for a specific project - working/input/output
    directories, path to scripts and other configs, sample names,
    options for specific rules in the pipeline, etc.
    locus.cluster_config.yaml: cluster configuration file for snakemake
    --cluster-config option. specifically for NIAID Locus HPC which uses
    UGE. (sets parameters for qsub command for each rule's job, and
    which environment modules to use)

-   locus_submit_vp.sh: batch job submit script for running the pipeline
    on Locus scripts: see scripts README

## Packages and versions (also found in locus_config.yaml) for all rules in pipeline

#### rule genomad: <https://portal.nersc.gov/genomad/installation.html>

-   genomad/1.5.2 -genomad database (genomaddb):
    /hpcdata/bio_data/genomad/genomad_db/ version 1.3 is what we are
    currently using. There is a updated version (1.7). The databases can
    be download directly from <https://zenodo.org/records/8339387> or
    with the command <genomad download-database> below

#### rule verse_genomad:<https://github.com/qinzhu/VERSE>

-   Verse/0.1.5 -Python/3.9.5-GCCcore-10.3.0

#### rule checkV: <https://bitbucket.org/berkeleylab/checkv/src/master/>

-   checkv/0.8.1 -Python/2.7.15-foss-2018b -checkv database (checkvdb):
    /hpcdata/bio_data/checkv/checkv-db-v1.5

#### rule checkv_filter: <https://github.com/lh3/seqtk>

-   Seqtk/1.3.r106

#### rule bbtools_dedupe:

<https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/installation-guide/>

-   bbmap/38.90

#### rule mmseqs: <https://github.com/soedinglab/MMseqs2> -mmseqs2/14.7e284

-   Python/3.9.5-GCCcore-10.3.0

#### rule vOTU

-   Python/3.9.5-GCCcore-10.3.0

#### rule vs4dramv: <https://github.com/jiarong/VirSorter2>

-   virsorter/2.2.3-Python-3.8.10

#### rule dramv and rule amgs: <https://github.com/WrightonLabCSU/DRAM>

-   dram/1.4.6 -dram database:
    /hpcdata/bcbb/shared/microbiome_share/modules/dram/1.4.6_ms

#### rule verse_dramv

-   Python/3.9.5-GCCcore-10.3.0 -liftoff/1.6.3-Python-3.9.12

#### rule iphop: <https://bitbucket.org/srouxjgi/iphop/src/main/>

-   iphop/1.3.0 -iphop database:
    /hpcdata/bio_data/iphop/Sept_2021_pub_rw/

#### rule diamond: <https://github.com/bbuchfink/diamond>

-   diamond/2.0.15.153 -diamond database (diamonddb)v0.9.28.129:
    /hpcdata/bio_data/DIAMOND_ncbi_db/nr.diamond.dmnd

#### rule gene_tables

-   Python/3.9.5-GCCcore-10.3.0
