# virome-pipeline

written by Lauren Krausfeldt & Poorani Subramanian - bioinformatics@niaid.nih.gov

This is a pipeline for exploring DNA viruses and bacteriophages in metagenomic samples. (still WIP)

### Files

- [_Snakefile_](Snakefile): pipeline script (reads in configs, commands for each pipeline step/rule)
  - [_cluster_setup.smk_](cluster_setup.smk): helper script for reading in cluster config file
- [_project_config.yaml_](project_config.yaml): for snakemake `--configfile` option.  config file with details for a specific project - working/input/output directories, path to scripts and other configs, sample names, options for specific rules in the pipeline, etc.
- [_locus.cluster_config.yaml_](locus.cluster_config.yaml): cluster configuration file for snakemake `--cluster-config` option. specifically for NIAID Locus HPC which uses UGE. (sets parameters for `qsub` command for each rule's job, and which environment modules to use)
- [_locus_submit_vp.sh_](locus_submit_vp.sh): batch job submit script for running the pipeline on Locus

## Running the Pipeline

### Inputs

The inputs to the pipeline are assembled contigs/scaffolds - one fasta file per sample.  They should be located in (or symlinked to) a single directory, and the filenames should start with a unique per-sample name.

### To run on Locus

1. Clone this repo locally (may need to give username and password):
```
git clone https://github.niaid.nih.gov/bcbb/virome-pipeline
```

2. Copy over the project config file _project_config.yaml_ and submit script _locus_submit_vp.sh_ to your project working directory, and edit both with the details for your specific project.  
   - for the submit script, the main items to edit are:
     1. path to the project config file
     2. the arguments for the `snakemake` command at the bottom of the script (see comments in the script)

3. Submit the job script:
  ```
  qsub ./locus_submit_vp.sh
  ```
  
4. Success?

## Notes

- This is tested to run on Locus.  However, it would be easy to adapt to another HPC that uses environment modules by making your own cluster config file (with the correct module names and job parameters), and your own job submit script (in particular modifying the `$clustercmd` for whatever job scheduler your HPC uses).
- In the future, we will work on making it more general (perhaps using conda or a containerized workflow instead of environment modules)
- Also, adding additional steps for specialized analysis and making the pipeline more flexible.
