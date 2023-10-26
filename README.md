# virome-pipeline

written by Lauren Krausfeldt & Poorani Subramanian - bioinformatics@niaid.nih.gov

## Description

This is a pipeline for exploring viruses (ssDNA, dsDNA phage, and giant DNA viruses) and viral diversity in metagenomes. The pipeline accepts assembly (.fasta) and binary alignment map (.bam) files as input. These files are produced from the WGSA2 pipeline in Nephele<sup>1</sup>. The output of this pipeline provides viral genomes found in the metagenome assembly, their taxonomy and level of completeness, viral functional genes and their abundances, and vOTU abundances and their host taxonomy. 

The pipeline first searchs for viral genomes using Genomad<sup>2</sup>, which also provides viral taxonomy and functional classification of each viral genomes. The viral genomes are also functionally classified with dramv<sup>3</sup> and (optionally) diamond<sup>4</sup> using the nr database. Gene abundances per sample are produced from these outputs using verse<sup>5</sup>. From here, the user has the option to filter the resulting sequences based on completeness using checkv<sup>6</sup>. Either the output of genomad or checkv is used to cluster viral genomes with bbtools dedupe<sup>7</sup> and mmseqs<sup>8</sup> to produce vOTUs<sup>9</sup>. Finally, abundances and host taxonomy of vOTUs are produced. 

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

## References

1. https://www.protocols.io/view/wgsa2-workflow-a-tutorial-n92ldm98xl5b/v1
2. Camargo, A. P., Roux, S., Schulz, F., Babinski, M., Xu, Y., Hu, B., ... & Kyrpides, N. C. (2023). Identification of mobile genetic elements with geNomad. Nature Biotechnology, 1-10. doi: [10.1038/s41587-023-01953-y](https://doi.org/10.1038/s41587-023-01953-y).   
3. Shaffer, M., Borton, M. A., McGivern, B. B., Zayed, A. A., La Rosa, S. L., Solden, L. M., ... & Wrighton, K. C. (2020). DRAM for distilling microbial metabolism to automate the curation of microbiome function. Nucleic acids research, 48(16), 8883-8900. doi: [10.1093/nar/gkaa621](https://doi.org/10.1093/nar/gkaa621).
4. Buchfink, B., Reuter, K., & Drost, H. G. (2021). Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature methods, 18(4), 366-368. doi: [10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x).
5. Zhu, Q., Fisher, S. A., Shallcross, J., & Kim, J. (2016). VERSE: a versatile and efficient RNA-Seq read counting tool. bioRxiv, 053306.  doi: [10.1101/053306](https://doi.org/10.1101/053306).
6. Nayfach, S., Camargo, A. P., Schulz, F., Eloe-Fadrosh, E., Roux, S., & Kyrpides, N. C. (2021). CheckV assesses the quality and completeness of metagenome-assembled viral genomes. Nature biotechnology, 39(5), 578-585. doi: [10.1038/s41587-020-00774-7](https://doi.org/10.1038/s41587-020-00774-7).
7. https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/dedupe-guide/
8. Steinegger, M., & Söding, J. (2017). MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature biotechnology, 35(11), 1026-1028. doi: [10.1038/nbt.3988](https://doi.org/10.1038/nbt.3988).
9. Roux, S., Adriaenssens, E. M., Dutilh, B. E., Koonin, E. V., Kropinski, A. M., Krupovic, M., ... & Eloe-Fadrosh, E. A. (2019). Minimum information about an uncultivated virus genome (MIUViG). Nature biotechnology, 37(1), 29-37. doi: [10.1038/nbt.4306](https://doi.org/10.1038/nbt.4306).
