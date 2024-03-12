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
    
## How it works
   
    
### Inputs
    
Input to the whole pipeline are an alignment file (.bam) and associated assemblies (.fasta).  
    
### Setup
    
Currently, the set up for the snakemake is that the input folder contains all .bam and .fasta file. The outputs all go into the virome_output folder. 

In the virome_output folder, a new folder gets made for each sample. 

For each sample, a folder gets made for each rule that is run on individual samples. 

Rules genomad, checkv, dramv, vs2, amgs, diamond, and verse are all run separately on each sample, so there will be separate folder for each of these results in each sample folder. 

Rules Bbtools_dudupe, mmseqs, vOTU run on combined data from all of the samples so those folders will be made in the virome_output folder. 

### Steps in the pipeline

### *Viral discovery:*

#### Rule genomad: (always run)

The tool Genomad is used to search metagenomic assemblies for viral genes and determine if a contig in an assembly is a full or partial genome. 

*Input:* fasta file (metagenome assembly) 

*User options:*
- Default (default)
- Conservative
- Relaxed

*Outputs:* viral genomes in fasta file 

#### Rule checkV (always run)
    
This rule cannot run before genomad is complete. 
CheckV is a tool that provide “quality” information about viral genomes found from gemonad.

*Input:* The output of genomad, which is a fasta file of viral sequences
*Output:* A tsv file showing quality metrics for each sequence. 

#### Rule checkV_filter (always run) 

This rule cannot run until checkV is complete. 
This allows the user to add more filtering of their viral contigs before any more processing of the data. 

*Input:* the output of genomad, which is a .fasta file of viral sequences. 
*User options:*
- All: Undetermined, low, medium, high, and complete (default)
- Medium quality, High quality, and complete 
- High quality, and complete 
- Complete 

*Output:*
Fasta file containing viral genomes that have been filtered for genome completeness. 

#### Rule verse_genomad (always run)

Verse cannot complete until genomad is finished. 
    
This rule extracts read counts from the bam files to get abundances of genes and all viral genomes from genomad’s output. 

*Input:* 
- original bam file that coincides with the metagenomic assembly
- output of genomad: gene information
- output of genomad: virus fasta file
- output of genomad: viral summary
- output of genomad: header map 

*Output:* Viral abundances and gene abundances from genomad output for each sample.

#### Rule dramV (always run)
    
This cannot run until checkv_filter finished.

*Input:* viral genomes from genomad output (fasta file) 

*Output:* 
- .tsv file with annotations from genes based on dramV 
- .tsv file with abundances of genes called by chosen database 
    
### *vOTU generation*
    
The following rules require all samples to have finished running through genomad and checkv_filter. 

#### Rule bbtools_dedupe (always run) 
    
Genomad and checkv_filter must complete before this starts. 
This rule collects all of the virus sequences from all samples submitted and combines them into one fasta file. Then the viruses that are 100% the same will be clustered. The longest sequence is retained. 

*Input:* Viral fasta sequences from all samples submitted into one fasta file. 
*Output:* Deduplicate viral fasta sequences from all samples in one fasta file.

#### Rule mmseqs (always run) 

Bbtools_dedupe rules must complete before this starts.
This rule clusters the deduplicated viral sequences from all samples to produce representative vOTU sequences. 

*Input:* deduplicated sequences from bbtools_dedupe (fasta file) 
*Output:* fasta file of vOTUs 

#### Rule vOTU (always run)
    
Cannot run until mmseqs and verse are finished. 
This rule combines the outputs from rule verse and provided abundances for vOTUs. It also merged the genomad summary information (containing taxonomy info) with the resulting vOTUs. 

*Input:* 
- tsv file that is an output from mmseqs
- virus counts from verse 
- .tsv genomad output summary information

*Output:*
- Abundance matrix of cpm values for each vOTU
- .tsv vOTU summary file 

#### Rule iphop (Optional)
    
Mmseqs must complete before this starts. 
This rule predicts host of virus. 

Input: .fasta output from mmseqs 

*Output:* .csv file with host information about all vOTUs.
    
### *Additional functional annotation (Whole step optional)*

#### Rule diamond: Optional
    
This cannot run until rule checkv_filter has finished. 

*Input:* protein sequences (fasta) from genomad output
*Output:* .tsv file with each genes and top hit in diamond database 

#### Rule vs24dramv and rule amgs:  Optional
    
This cannot run until checkv_filter has finished. Rule amgs will not run unless vs2dramv is selected and if it is, it must also finish. 

*Input:* viral genomes from genomad output (fasta file) 
*Output:*
- .tsv file with annotations from genes based on dramV 
- .tsv file of amg abundances (CPM) across samples


  
