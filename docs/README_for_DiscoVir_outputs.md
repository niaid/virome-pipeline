# README of DiscoVir Outputs

Thanks for using [Nephele's](https://nephele.niaid.nih.gov) DiscoVir virome pipeline!  See our [pipeline details page](https://nephele.niaid.nih.gov/pipeline_details/discovir) for more information, and for additional questions and help, please email us at nephelesupport@nih.gov.

There are two main folder types in the output of the DiscoVir pipeline:

- Sample folders: per-sample output
- Combined folders: outputs of analysis on sequences combined from across the entire dataset & tables merged from the individual sample results.

**Combined folders**

**Sample folders**

Each sample has an individual folder which contains subfolders for each step in the pipeline that analyzes the individual sample's sequences.  Here we will highlight the final output of each step as well as other important files.  For a complete list of the outputs of each tool (including any intermediate files), please see the individual tool's documentation linked below.  

The subfolders inside the sample folder are:

- ***genomad**:* outputs of [geNomad](https://portal.nersc.gov/genomad/).

  - *{sample}.genomad.log*: log file to check for all (STDOUT and STDERR) messages from this pipeline step.
  - *{sample}_summary*: final output of genomad
    - *{sample}\_summary/{sample}\_virus_summary.tsv*: summary table of viral sequences identified
    - *{sample}\_summary/{sample}\_virus.fna*: FASTA of viral sequences (with host trimmed, if necessary)
    - *{sample}\_summary/{sample}\_virus_genes.tsv*: summary table of genes predicted from viral sequences
    - *{sample}\_summary/{sample}\_virus_proteins.faa*: amino acid/protein FASTA of viral genes
  - *abund\_genomad/{sample}\_virus.count.CDS.cpm.txt:* abundance estimates of the viral sequences.
    - readcounts are estimated using the ht-seq algorithm from the tool [VERSE](https://anaconda.org/bioconda/verse) , along with estimated copies per million (cpm), and reads per kilobase (rpk)
    - cpm is calculated the same way as tpm, but here we use as a generic term in case of DNA-seq.
  - *abund\_genomad/{sample}\_virus_genes.count.CDS.cpm.txt:* abundance estimates of viral genes

- **checkv:** outputs of [CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/)
  - *{sample}.checkv.log:* log file
  - *quality_summary.tsv:* summary of quality of all viral sequences
  - *combined.fna:* viral sequences identified by CheckV (with headers given by CheckV)
  - *checkv_filtered_genomad_viruses.fna:* viral sequences filtered based on the `checkv quality` user option.  this file uses the original genomad headers which allows the pipeline to more easily track the sequences through analysis.

- **dramv:** outputs of running [DRAM-v](https://github.com/WrightonLabCSU/DRAM/wiki/4b.-Interpreting-the-Results-of-DRAM-v) annotate step only using the viral sequences predicted by geNomad and optionally filtered by CheckV ( *checkv/checkv_filtered_genomad_viruses.fna*).  DRAM-v finds genes and annotates them.

  - *{sample}.dramv.log*: log file
  - *dramv-annotate/annotations.tsv:* table of DRAM-v gene annotations
  - *dramv-annotate/genes.{faa,fna,gff}*: AA and nucleotide FASTA files as well as associated gff file for the genes.  the coordinates in the gff are vis-a-vis the original geNomad sequences found in *checkv/checkv_filtered_genomad_viruses.fna*.

  - *abund_dramv/{sample}\_dramv.count.gene.cpm.txt:* abundance estimates of genes found by DRAM-V

- **amgs (optional):** if the `run_amgs` user option is chosen, this folder has the [DRAM-v](https://github.com/WrightonLabCSU/DRAM/wiki/4b.-Interpreting-the-Results-of-DRAM-v) predicted AMGs (auxiliary metabolic genes).  To predict the AMGs, the pipeline first runs [VirSorter2](https://github.com/jiarong/VirSorter2) on *checkv_filtered_genomad_viruses.fna* to generate the table needed by DRAM-v to predict AMGs.  VirSorter2 filters out many sequences identified as viral by other tools, so this is an optional step.

  - *vs2:* output of [VirSorter2](https://github.com/jiarong/VirSorter2)
    - *vs2/{sample}.vs2.log:* log file
    - *vs2/for-dramv/final-viral-combined-for-dramv.fna & vs2/for-dramv/viral-affi-contigs-for-dramv.tab:* files used by DRAM-v to predict AMGs
  - *{sample}.dramv_amgs.log:* DRAM-v log file
  - *dramv-annotate/:* directory of DRAM-v gene annotations - using all DRAM-v databases
  - *dramv-distill:* directory containing the output of DRAM-v's distill step which identifies AMGs
    - *dramv-distill/amg_summary.tsv:* table of potential AMGs 
  - *abund_amgs/{sample}_amgs.count.gene.cpm.txt:* abundance estimates of genes

- **diamond (optional):** output of annotating the geNomad-predicted genes by aligning sequences with diamond to NCBI's nr database.  

  - *CA1.nr.diamond.tsv:* table of top alignments for gene sequences with NCBI nr accession number.  For full explanation of all columns see the [NCBI BLAST format table under *outfmt*](https://www.ncbi.nlm.nih.gov/books/NBK279684/#appendices.Options_for_the_commandline_a) (DIAMOND uses the BLAST output format).