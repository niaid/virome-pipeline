from os.path import join as pjoin
import os
import snakemake.utils as sm

### cluster config
scriptdir = config["scriptdir"] ## virome-pipeline repo/script directory
include: pjoin(scriptdir, "cluster_setup.smk")

############ project config #################

## directories
workdir: config["workdir"]  ## working directory
OUT = config["outdir"] ## directory storing all output
SOUT = pjoin(OUT, "{sample}") ## subdirectory for each sample's output (including logs)
LOGS = pjoin(OUT, "logs") ## log directory for logfiles of rules involving more than one sample
IN = config["indir"] ## directory storing all assemblies/input
VOTU_DIR = pjoin(OUT, "vOTUs")
VOTU_CL_DIR = pjoin(VOTU_DIR, "vOTU_clustering")
HOST_DIR = pjoin(OUT, "vOTU_Host_prediction_iphop")


## inputs
SAMPLES = list(config["samples"].split())
if SAMPLES is None:
    sys.exit("Sample prefixes are needed")

################ RULES ###########################
onerror:
     print("An error occurred in the virome pipeline.")
     shell("mail -s 'Error in virome pipeline.' {config[email]} < {log}")

localrules: all, createsampledir


####### SETUP RULES #######

rule createsampledir:
    """setup output directories"""
    output: pjoin(SOUT, "logs/snakefake")
    params: outdir = SOUT
    shell:"""
    mkdir -p {params.outdir}/logs
    touch {output}
    """


########## PIPELINE RULES ########

genomad_filter = { "default": "",
                   "relaxed": "--relaxed",
                   "conservative": "--conservative"
                   }
rule genomad:
    threads: clust_conf["genomad"]["threads"]
    envmodules: clust_conf["genomad"]["modules"]
    input: fake = ancient(rules.createsampledir.output),
           assembly = pjoin(IN, "{sample}" + config["assembly_suffix"])
    params: outdir = pjoin(SOUT, "genomad"),
            renamed_assembly = pjoin(SOUT, "genomad", "{sample}" + ".fasta"),
            genomad_filter = genomad_filter[config["genomad_filter"]]
    log:    pjoin(SOUT, "genomad", "{sample}" + ".genomad.log")
    output: fna = pjoin(SOUT, "genomad", "{sample}_summary", "{sample}_virus.fna"),
            genes = pjoin(SOUT, "genomad", "{sample}_summary", "{sample}_virus_genes.tsv"),
            proteins =  pjoin(SOUT, "genomad", "{sample}_summary", "{sample}_virus_proteins.faa"),
            summary = pjoin(SOUT, "genomad", "{sample}_summary","{sample}_virus_summary.tsv"),
            assembly_headermap = pjoin(SOUT, "genomad", "{sample}" + "_assembly_headermap.txt")
    shell:"""
    ## cleanup possible previous run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    ## genomad search for viral contigs
    genomad --version
    genomad --version 1>>{log}



    ## rename assembly to sample.fasta so genomad uses sample as output file prefix
    ## also rename headers to include sample name, so we can track contigs when we combine all for clustering
    renamed_assembly=$(realpath {params.renamed_assembly})

    sed '/^>/ s!\([^ ]*\)!\\1_{wildcards.sample}!' < {input.assembly} >$renamed_assembly

    grep -e "^>" {input.assembly} | sed -r 's/^.{{1}}//' >{params.outdir}/oldheaders.txt
    grep -e "^>" $renamed_assembly | sed -r 's/^.{{1}}//' >{params.outdir}/newheaders.txt
    paste {params.outdir}/oldheaders.txt {params.outdir}/newheaders.txt >{output.assembly_headermap}
    rm {params.outdir}/oldheaders.txt
    rm {params.outdir}/newheaders.txt


    ## run genomad

    echo "Looking for viruses in {wildcards.sample} using geNomad.  See log file {log}."

    genomad end-to-end {params.genomad_filter} --enable-score-calibration --cleanup --threads {threads} $renamed_assembly {params.outdir} {config[genomaddb]} 1>>{log}

    echo "geNomad finished for sample {wildcards.sample}"

    rm $renamed_assembly

    """

rule abund_genomad:
    threads: clust_conf["abund_genomad"]["threads"]
    envmodules: *clust_conf["abund_genomad"]["modules"]
    input:  genes = rules.genomad.output.genes,
            fna = rules.genomad.output.fna,
            virus = rules.genomad.output.summary,
            headermap = rules.genomad.output.assembly_headermap,
            bam = ancient(pjoin(IN, "{sample}" + config["bam_suffix"]))
    params: outdir = pjoin(SOUT, "genomad", "abund_genomad"),
            prefix_genes = pjoin(SOUT, "genomad", "abund_genomad", "{sample}_virus_genes.count"),
            counts_only_genes = pjoin(SOUT, "genomad", "abund_genomad", "{sample}_virus_genes.count.CDS.txt"),
            prefix_virus = pjoin(SOUT, "genomad", "abund_genomad", "{sample}_virus.count"),
            counts_only_virus = pjoin(SOUT, "genomad", "abund_genomad", "{sample}_virus.count.CDS.txt")
    log:    pjoin(SOUT, "genomad", "abund_genomad", "{sample}" + ".abund_genomad.log")
    output: readcounts_genes = pjoin(SOUT, "genomad", "abund_genomad", "{sample}_virus_genes.count.CDS.cpm.txt"),
            readcounts_virus = pjoin(SOUT, "genomad", "abund_genomad", "{sample}_virus.count.CDS.cpm.txt")
    shell:"""
    ## cleanup possible previous run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    ## estimate gene abundances with verse
    verse -v
    verse -v >>{log}

    echo "Estimating abundances of genes from geNomad for {wildcards.sample} with verse. See log file {log}." 

    ## make gff from virus genes files
    python3 {config[scriptdir]}/scripts/genomad_genes2gff.py {input.genes} -m {input.headermap} >{params.prefix_genes}.gff
    verse -a {params.prefix_genes}.gff -o {params.prefix_genes} -g ID -z 1 -t CDS -l -T {threads} {input.bam} 1>>{log}

    python3 {config[scriptdir]}/scripts/calc_cpm.py {params.counts_only_genes} >{output.readcounts_genes}

    echo "Estimating viral abundances for {wildcards.sample} with verse.  See log file {log}." 

    ## make gff from virus summary; use default -z 1 instead of -z 5
    python3 {config[scriptdir]}/scripts/genomad_virus2gff.py {input.virus} -m {input.headermap} >{params.prefix_virus}.gff
    verse -a {params.prefix_virus}.gff -o {params.prefix_virus} -g ID -z 1 -t CDS -l -T {threads} {input.bam} 1>>{log}

    python3 {config[scriptdir]}/scripts/calc_cpm.py {params.counts_only_virus} >{output.readcounts_virus}


    """


rule checkv:
    threads: clust_conf["checkv"]["threads"]
    envmodules: *clust_conf["checkv"]["modules"]
    input: rules.genomad.output.fna
    params: outdir = pjoin(SOUT, "checkv")
    log: pjoin(SOUT, "checkv", "{sample}" + ".checkv.log")
    output: fna = pjoin(SOUT, "checkv", "combined.fna"),
            qsum = pjoin(SOUT, "checkv", "quality_summary.tsv")
    shell:"""
    ## run checkv to qc genomad results and trim host regions left at the end of proviruses

    ## cleanup possible previous run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    checkv | head -n1
    checkv | head -n1 >>{log}

    echo "Checking quality and completeness of viral genomes in {wildcards.sample} with checkV.  See log file for output and any errors {log}."

    checkv end_to_end {input} {params.outdir} -d {config[checkvdb]} -t {threads} --quiet 1>>{log} 2>>{log}

    cat {params.outdir}/proviruses.fna {params.outdir}/viruses.fna >{output.fna}

    echo "CheckV finished for {wildcards.sample}."

    """

def checkv_q_validate(checkv_q):
    """return grep string for filtering checkv contigs"""
    if checkv_q not in ['all', 'medium', 'high', 'complete']:
        raise ValueError(f"checkv_q value in config is '{checkv_q}'.  Should be one of 'all', 'medium', 'high', 'complete'. Exiting.")
    return(checkv_q)


rule checkv_filter:
    threads: clust_conf["checkv_filter"]["threads"]
    envmodules: clust_conf["checkv_filter"]["modules"]
    input: fna = rules.genomad.output.fna,
           qsum = rules.checkv.output.qsum
    params: outdir = pjoin(SOUT, "checkv"),
            checkv_q = checkv_q_validate(config["checkv_q"]),
            tempclist = pjoin(SOUT, "checkv", "tempclist.txt"),
            quality = pjoin(SOUT, "checkv", "checkv_quality_genomad_viruses.fna")
    output: pjoin(SOUT, "checkv", "checkv_filtered_genomad_viruses.fna")
    shell:"""
    ## use checkv quality assessment to filter viral contigs for clustering

    ## cleanup possible previous run
    rm -rf {output}

    echo "Filtering viral genomes by quality for {wildcards.sample}."

    if [ "{params.checkv_q}" = "complete" ]; then
       grep -e "Complete" {input.qsum} | awk '{{ print $1 }}' >{params.tempclist}
    elif [ "{params.checkv_q}" = "high" ]; then
       grep -e "Complete" -e "High-quality" {input.qsum} | awk '{{ print $1 }}' >{params.tempclist}
    elif [ "{params.checkv_q}" = "medium" ]; then
       grep -e "Complete" -e "High-quality" -e "Medium-quality" {input.qsum} | awk '{{ print $1 }}' >{params.tempclist}
    fi

    if [ "{params.checkv_q}" != "all" ]; then
        seqtk subseq {input.fna} {params.tempclist} > {params.outdir}/temp && mv {params.outdir}/temp {params.quality}

        rm -fv {params.tempclist}
    else
        cp {input.fna} {params.quality}
    fi
    
     seqtk seq -L {config[vs_min_length]} {params.quality} > {output}


    """



rule bbtools_dedupe:
    threads: clust_conf["bbtools_dedupe"]["threads"]
    envmodules: clust_conf["bbtools_dedupe"]["modules"]
    input: expand(rules.checkv_filter.output, sample=SAMPLES)
    params: outdir = pjoin(VOTU_CL_DIR, "bbtools_dedupe"),
            input_ctgs = pjoin(VOTU_CL_DIR, "bbtools_dedupe", "all_input_contigs.fasta"),
	    log = pjoin(VOTU_CL_DIR, "bbtools_dedupe", "log.txt"),
	    stats = pjoin(VOTU_CL_DIR, "bbtools_dedupe", "cluster_stats.txt")
    log:  pjoin(VOTU_CL_DIR, "bbtools_dedupe", "bbtools_dedupe.log")
    output: unique_seqs = pjoin(VOTU_CL_DIR, "bbtools_dedupe", "unique_seqs.fasta")
    shell:"""
    ## run bbtools_dedupe on all contigs to remove duplicates

    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    dedupe.sh --version

    echo "Starting vOTU generation. See log file {log}."
  
    echo "Collecting viral genomes from all samples." 

    cat {input} >{params.input_ctgs}

    echo "Deduplicating viral genomes."

    dedupe.sh in={params.input_ctgs} out={output.unique_seqs} csf={params.stats} minscaf={config[vs_min_length]} \
	mergenames=t ex=f usejni=t threads={threads} 1>>{log}


    """

rule mmseqs:
    threads: clust_conf["mmseqs"]["threads"]
    envmodules: *clust_conf["mmseqs"]["modules"]
    input: rules.bbtools_dedupe.output.unique_seqs
    params: outdir = pjoin(VOTU_CL_DIR, "mmseqs"),
            DB_dir = pjoin(VOTU_CL_DIR, "mmseqs", "DB"),
            DB = pjoin(VOTU_CL_DIR, "mmseqs", "DB/DB"),
	    DB_clu = pjoin(VOTU_CL_DIR, "mmseqs", "DB/DB_clu"),
	    DB_clu_seq = pjoin(VOTU_CL_DIR, "mmseqs", "DB/DB_clu_seq"),
	    DB_clu_fasta = pjoin(VOTU_CL_DIR, "mmseqs", "cluster_seqs.fasta"),
	    DB_clu_rep = pjoin(VOTU_CL_DIR, "mmseqs", "DB/DB_clu_rep")
    log: pjoin(VOTU_CL_DIR, "mmseqs", "mmseqs2.log")
    output: DB_clu_rep_fasta = pjoin(VOTU_CL_DIR, "mmseqs", "representative_seqs.fasta"),
            DB_clu_tsv = pjoin(VOTU_CL_DIR, "mmseqs", "DB_clu.tsv"),
            flat_DB_clu_tsv = pjoin(VOTU_CL_DIR, "mmseqs", "flat_DB_clu.tsv"),
            renamed_DB_clu_rep_fasta = pjoin(VOTU_CL_DIR, "mmseqs", "representative_seqs.renamed.fasta"),
            vOTU_fasta = pjoin(VOTU_DIR, "vOTU_sequences.fasta")
    shell:"""
    ## run mmseqs on all deduped genomes

    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    mmseqs version

    mkdir {params.DB_dir}

    echo "Clustering unique viral genomes with 95% identity and 85% coverage to generate vOTUs.  See log file {log}."

    mmseqs createdb {input} {params.DB} 1>>{log}

    mmseqs cluster {params.DB} {params.DB_clu} $TMPDIR --cov-mode 1 -c 0.85 \
	--min-seq-id 0.95 --cluster-mode 2 --threads {threads} 1>>{log}

    mmseqs createtsv {params.DB} {params.DB} {params.DB_clu} {output.DB_clu_tsv} \
        --threads {threads} 1>>{log}

    mmseqs createseqfiledb {params.DB} {params.DB_clu} {params.DB_clu_seq} \
        --threads {threads} 1>>{log}

    mmseqs result2flat {params.DB} {params.DB} {params.DB_clu_seq} \
        {params.DB_clu_fasta} 1>>{log}

    mmseqs createsubdb {params.DB_clu} {params.DB} {params.DB_clu_rep} 1>>{log}

    mmseqs convert2fasta {params.DB_clu_rep} {output.DB_clu_rep_fasta} 1>>{log}

    ## rename repseq to only use first dup sequence from bbtools_dedupe
    # flatten clu.tsv file so each row is one repseq and one seq
    python3 {config[scriptdir]}/scripts/flatten_mmseqs_tsv.py {output.DB_clu_tsv} | sort -u > {output.flat_DB_clu_tsv}
    # rename representative_seqs.fasta
    sed '/^>/ s!>[^>]*!!2g' {output.DB_clu_rep_fasta} >{output.renamed_DB_clu_rep_fasta}

    cp {output.renamed_DB_clu_rep_fasta} {output.vOTU_fasta}

    echo "Generation of vOTUs finished."


    """

rule votu:
    threads: clust_conf["votu"]["threads"]
    envmodules: *clust_conf["votu"]["modules"]
    input: mmseqs = rules.mmseqs.output.flat_DB_clu_tsv,
           abund = expand(rules.abund_genomad.output.readcounts_virus, sample=SAMPLES),
           summ = expand(rules.genomad.output.summary, sample=SAMPLES)
    params: outdir = pjoin(VOTU_DIR),
            filelist = pjoin(VOTU_DIR, "abundfiles.txt"),
            summlist = pjoin(VOTU_DIR, "summaryfiles.txt"),
            pythonpath = clust_conf["votu"]["pythonpath"],
            mapping = config["mapping_file"],
            obs_met = pjoin(VOTU_DIR, "obs_met.tsv"),
            tempbiom = pjoin(VOTU_DIR, "temp.biom")
    output: votu = pjoin(VOTU_DIR, "vOTU_table_cpm.tsv"),
            gmdanno = pjoin(VOTU_DIR, "vOTUs_genomad_virus_summary.tsv"),
            biom = pjoin(VOTU_DIR, "vOTU_cpm.biom")

    shell:"""

    ## make vOTU table
    python3 --version

    ## cleanup possible previous failed run
    rm -f {output.votu} {output.gmdanno} {params.filelist} {params.summlist} {params.obs_met} {params.tempbiom}

    echo "Creating vOTU abundance table." 

    echo "{input.abund}" >{params.filelist}
    tr " " "\\n" <{params.filelist} >{params.outdir}/temp && mv {params.outdir}/temp {params.filelist}
  
    python3 {config[scriptdir]}/scripts/make_votu_table.py -v cpm {input.mmseqs} {params.filelist} >{output.votu}

    ## collate genomad annotations
    echo "{input.summ}" | tr " " "\\n" >{params.summlist}
    python3 {config[scriptdir]}/scripts/repseq_genomad_annotations.py {output.votu} \
        {params.summlist} >{output.gmdanno}

    # make krona charts
    mkdir -p {params.outdir}/temp1
    python3 {config[scriptdir]}/scripts/format_genomad_and_votu_for_krona.py {output.votu} {output.gmdanno} {params.outdir}/temp1
    ktImportText -o {params.outdir}/vOTU.krona.html {params.outdir}/temp1/*.txt

    ## make biom file
    sed '0,/seq_name/{{s/seq_name/\#repseq/}}' {output.gmdanno} | awk -F $'\t' '{{ print $1"\t"$11 }}' >{params.obs_met}
    export PYTHONPATH={params.pythonpath}
    {params.pythonpath}/bin/biom convert -i {output.votu} -o {params.tempbiom} --to-json --table-type "OTU table"
    {params.pythonpath}/bin/biom add-metadata -i {params.tempbiom} -o {output.biom} \
         --observation-metadata-fp {params.obs_met} --sc-separated taxonomy --output-as-json
    if [ "{params.mapping}" != "None" ]; then
        {params.pythonpath}/bin/biom add-metadata -i {output.biom} -o {params.tempbiom} \
         --sample-metadata-fp {params.mapping} --output-as-json
    fi

  
    ## cleanup
    rm -f {params.filelist} {params.summlist} {params.obs_met} {params.tempbiom}
    rm -rf {params.outdir}/temp1

    """

rule vs4dramv:
    threads: clust_conf["vs4dramv"]["threads"]
    envmodules: clust_conf["vs4dramv"]["modules"]
    input: rules.checkv_filter.output
    params: outdir = pjoin(SOUT, "amgs", "vs2")
    log: pjoin(SOUT, "amgs/vs2", "{sample}" + ".vs2.log")
    output: fasta = pjoin(SOUT, "amgs/vs2/for-dramv/final-viral-combined-for-dramv.fa"),
            tab = pjoin(SOUT, "amgs/vs2/for-dramv/viral-affi-contigs-for-dramv.tab")
    shell:"""
    ## run virsorter to produce input for DRAM-v
    virsorter --version

    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    echo "Running all viral genomes from {wildcards.sample} through Virsorter2.0 to look for AMGs.  See log file {log}."

    virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off \
        --prep-for-dramv -i {input} -w {params.outdir} --include-groups dsDNAphage,ssDNA,NCLDV \
        --min-length {config[vs_min_length]} --min-score 0.5 --rm-tmpdir -j {threads} all 1>>{log} 2>>{log}

    echo "Virsorter2.0 finished for {wildcards.sample}."

    """

rule amgs:
    threads: clust_conf["amgs"]["threads"]
    envmodules: clust_conf["amgs"]["modules"]
    input: fasta = rules.vs4dramv.output.fasta,
           tab = rules.vs4dramv.output.tab
    params: annot_outdir = pjoin(SOUT, "amgs", "dramv-annotate"),
            distill_outdir = pjoin(SOUT, "amgs", "dramv-distill")
    log: pjoin(SOUT, "amgs", "{sample}" + ".dramv_amgs.log")
    output: genes = pjoin(SOUT, "amgs/dramv-annotate/genes.faa"),
            gff = pjoin(SOUT, "amgs/dramv-annotate/genes.gff"),
            summary = pjoin(SOUT, "amgs/dramv-distill/amg_summary.tsv"),
            scaffolds = pjoin(SOUT, "amgs", "dramv-annotate", "scaffolds.fna")

    shell:"""
    ## DRAM-v annotation of viral sequences for AMGS
    DRAM-setup.py version
    mmseqs -h | head -n 7

    ## cleanup possible previous failed run
    rm -rf {params.annot_outdir}
    rm -rf {params.distill_outdir}
    rm -f {log}

   echo "Annotating all viral genomes from {wildcards.sample} with dramv to look for AMGs. See log file {log}."
    
        DRAM-v.py annotate -i {input.fasta} \
        -v {input.tab} \
        -o {params.annot_outdir} --skip_trnascan \
        --threads {threads} --min_contig_size {config[vs_min_length]} 1>>{log} 2>>{log}
       
        DRAM-v.py distill -i {params.annot_outdir}/annotations.tsv \
        -o {params.distill_outdir} --log_file_path {log} 1>>{log} 2>>{log}

    echo "dramv for amgs from {wildcards.sample} finished."

    """

rule abund_amgs:
    threads: clust_conf["abund_amgs"]["threads"]
    envmodules: *clust_conf["abund_amgs"]["modules"]
    input:  gff = rules.amgs.output.gff,
            target = pjoin(IN, "{sample}" + config["assembly_suffix"]),
            reference = pjoin(SOUT, "amgs", "dramv-annotate", "scaffolds.fna"),
            bam = ancient(pjoin(IN, "{sample}" + config["bam_suffix"]))
    log:    pjoin(SOUT, "amgs", "abund_amgs", "{sample}" + ".abund_amgs.log")
    params: outdir = pjoin(SOUT, "amgs", "abund_amgs"),
            intermdir = pjoin(SOUT, "amgs", "abund_amgs", "intermediate_files"),
            prefix_genes = pjoin(SOUT, "amgs", "abund_amgs", "{sample}_amgs.count"),
            counts_only_genes = pjoin(SOUT, "amgs", "abund_amgs", "{sample}_amgs.count.gene.txt"),
            liftoff_gff = pjoin(SOUT, "amgs", "abund_amgs", "{sample}.amgs_genes.liftoff.gff")
    output: readcounts_genes = pjoin(SOUT, "amgs", "abund_amgs", "{sample}_amgs.count.gene.cpm.txt")
    shell:"""
    ## lift gene features from trimmed virsorter2 contigs to original with liftoff and estimate AMG abundances with verse
    liftoff -V
    verse -v

    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.intermdir}

    echo "Calculating AMG abundances.  For liftoff output and errors, see {log}." 

    # make chrom file
    grep -v -w "##gff-version" {input.gff} | grep -v -e "^#" | awk '{{ print $1 }}' | sort -u >{params.intermdir}/1.txt
    cat {params.intermdir}/1.txt | sed "s/\(.*\)_{wildcards.sample}.*/\\1/"  >{params.intermdir}/2.txt
    paste -d "," {params.intermdir}/1.txt {params.intermdir}/2.txt | sort -u >{params.intermdir}/chroms.txt

    # change CDS to gene in dramv gff as liftoff works with gene features
    sed $'s/\tCDS\t/\tgene\t/' < {input.gff} >{params.intermdir}/temp.gff

    ## run liftoff
    liftoff {input.target} {input.reference} -g {params.intermdir}/temp.gff -o {params.liftoff_gff} \
        -u {params.outdir}/unmapped_features.txt \
        -dir {params.outdir}/intermediate_files \
        -p 1 -chroms {params.intermdir}/chroms.txt \
        -exclude_partial -a 0.9 2>>{log} 1>>{log}

    ## remove double quotes
    sed 's/\"//g' {params.liftoff_gff} >{params.outdir}/temp.liftoff.gff && mv {params.outdir}/temp.liftoff.gff {params.liftoff_gff}

    ## run verse
    verse -a {params.liftoff_gff} -o {params.prefix_genes} -g ID -z 1 -t gene \
        -l -T 8 {input.bam} 1>>{log}

    ## calc cpms
    python3 {config[scriptdir]}/scripts/calc_cpm.py {params.counts_only_genes} >{output.readcounts_genes}
    rm {params.counts_only_genes}
    rm -rf {params.intermdir}


    """

rule dramv:
    threads: clust_conf["dramv"]["threads"]
    envmodules: clust_conf["dramv"]["modules"]
    input: fasta = rules.checkv_filter.output
    params: outdir = pjoin(SOUT, "dramv")
    log: pjoin(SOUT, "dramv", "{sample}" + ".dramv.log")
    output: genes = pjoin(SOUT, "dramv/dramv-annotate/genes.faa"),
            gff = pjoin(SOUT, "dramv/dramv-annotate/genes.gff"),
            scaffolds = pjoin(SOUT, "dramv/dramv-annotate/scaffolds.fna")

    shell:"""
    ## DRAM-v annotation of viral sequences
    DRAM-setup.py version
    mmseqs -h | head -n 7

    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    echo "Running dramv on all viral sequences from geNomad in {wildcards.sample} for functional annotation. See {log}." 
  
        DRAM-v.py annotate -i {input.fasta} \
        -o {params.outdir}/dramv-annotate --skip_trnascan \
        --threads {threads} --min_contig_size {config[vs_min_length]} 1>>{log} 2>>{log}

    echo "dramv for functional annotation on viral genomes from geNomad in {wildcards.sample} finished." 

    """

rule abund_dramv:
    threads: clust_conf["abund_dramv"]["threads"]
    envmodules: *clust_conf["abund_dramv"]["modules"]
    input:  gff = rules.dramv.output.gff,
            target = pjoin(IN, "{sample}" + config["assembly_suffix"]),
            reference =  rules.dramv.output.scaffolds,
            bam = ancient(pjoin(IN, "{sample}" + config["bam_suffix"]))
    params: outdir = pjoin(SOUT, "dramv", "abund_dramv"),
            intermdir = pjoin(SOUT, "dramv", "abund_dramv", "intermediate_files"),
            prefix_genes = pjoin(SOUT, "dramv", "abund_dramv", "{sample}_dramv.count"),
            counts_only_genes = pjoin(SOUT, "dramv", "abund_dramv", "{sample}_dramv.count.gene.txt"),
            liftoff_gff = pjoin(SOUT, "dramv", "abund_dramv", "{sample}.dramv_genes.liftoff.gff")
    log: pjoin(SOUT, "dramv", "abund_dramv", "{sample}" + ".abund_dramv.log")
    output: readcounts_genes = pjoin(SOUT, "dramv", "abund_dramv", "{sample}_dramv.count.gene.cpm.txt")
    shell:"""
    ## lift gene features from trimmed contigs to original with liftoff and estimate gene abundances with verse
    liftoff -V
    verse -v
    

    ## cleanup possible previous run
    rm -rf {params.outdir}
    mkdir -p {params.intermdir}

   echo "Calculating abundances of genes identified and annotated by DRAM-v using Liftoff and VERSE.  See log file for errors and output {log}." 

    # make chrom file
    grep -v -w "##gff-version" {input.gff} | grep -v -e "^#" | awk '{{ print $1 }}' >{params.intermdir}/1.txt
    cat {params.intermdir}/1.txt | sed 's/|.*//' | sed "s/\(.*\)_{wildcards.sample}.*/\\1/"  >{params.intermdir}/2.txt
    paste -d "," {params.intermdir}/1.txt {params.intermdir}/2.txt | sort -u >{params.intermdir}/chroms.txt
    

    ## rename CDS features to gene as liftoff only transfers genes. 
    sed $'s/\tCDS\t/\tgene\t/' < {input.gff} >{params.intermdir}/temp.gff

    ## reduce -p parallel to 1 as sometimes
    ## we get core dumps if the contigs are very large.
    ## redirect output to a log file otherwise it overwhelms the main log
    liftoff {input.target} {input.reference} -g {params.intermdir}/temp.gff -o {params.liftoff_gff} \
          -u {params.outdir}/unmapped_features.txt -dir {params.intermdir} -p 1 \
          -chroms {params.intermdir}/chroms.txt -exclude_partial -a 0.9 2>{log} \
          1>>{log}

    ## remove double quotes
    sed 's/\"//g' {params.liftoff_gff} >{params.outdir}/temp.liftoff.gff && mv {params.outdir}/temp.liftoff.gff {params.liftoff_gff}

    ## get reads counts
    verse -a {params.liftoff_gff} -o {params.prefix_genes} -g ID -z 1 -t gene -l -T {threads} {input.bam} 1>>{log}

    python3 {config[scriptdir]}/scripts/calc_cpm.py {params.counts_only_genes} >{output.readcounts_genes}
    rm {params.counts_only_genes}
    rm -rf {params.intermdir}

    """

rule gene_tables:
    threads: clust_conf["gene_tables"]["threads"]
    envmodules: clust_conf["gene_tables"]["modules"]
    input: abund = expand(rules.abund_dramv.output.readcounts_genes, sample=SAMPLES)
    params: outdir = pjoin(OUT, "gene_tables"),
            samp = SAMPLES,
            samplelist = pjoin(OUT, "gene_tables", "samplelist.txt"),
            workingdir = OUT,
            top_vogids = pjoin(OUT, "gene_tables", "top_dramv_vogdb_hits_cpm.tsv"),
            vogdb_heatmap = pjoin(OUT, "gene_tables", "dramv_vogdb_heatmap_cpm"),
            pythonpath = clust_conf["gene_tables"]["pythonpath"]
    output: pfam = pjoin(OUT, "gene_tables", "dramv_pfam_hits_cpm.tsv"),
            vogdb = pjoin(OUT, "gene_tables", "dramv_vogdb_hits_cpm.tsv"),
            kofam = pjoin(OUT, "gene_tables", "dramv_kofam_hits_cpm.tsv")
    shell:"""
    ## make abundance tables for dramv and diamond genes over all samples

    ## cleanup possible previous failed run
    rm -f {output.pfam} {output.vogdb} {output.kofam}
    mkdir -p {params.outdir}

    ## make input sample list for script
    echo "{params.samp}" >{params.samplelist}
    tr " " "\\n" <{params.samplelist} >{params.outdir}/temp && mv {params.outdir}/temp {params.samplelist}

    ## dramv-annotate abund tables

    echo "Collating abundances of VOGIDs, pfams, and kofams from dramv functional annotation." 

    python3 {config[scriptdir]}/scripts/dramv_genes_table.py {params.workingdir} {params.samplelist} -v cpm -c pfam_hits >{output.pfam}
    python3 {config[scriptdir]}/scripts/dramv_genes_table.py {params.workingdir} {params.samplelist} -v cpm -c ko_id -c kegg_hit >{output.kofam}
    python3 {config[scriptdir]}/scripts/dramv_genes_table.py {params.workingdir} {params.samplelist} -v cpm -c vogdb_id -c vogdb_hits >{output.vogdb}

    ## vogdb heatmap of top genes by prevalence and mean
    echo "Making heatmap of top VOG genes by prevalence and abundance." 
    python3 {config[scriptdir]}/scripts/top_genes.py {output.vogdb} -i "vogdb_hits" -m prev -n 25 >{params.top_vogids}
    export PYTHONPATH={params.pythonpath}
    python3 {config[scriptdir]}/scripts/plotnine_heatmap.py {params.top_vogids} {params.vogdb_heatmap} \
           -t "Top VOG genes" -d "vogdb_hits" -a "cpm" 
     

    ## rm {params.samplelist}

    """


rule amg_tables:
    threads: clust_conf["amg_tables"]["threads"]
    envmodules: clust_conf["amg_tables"]["modules"]
    input: amgs = expand(rules.abund_amgs.output.readcounts_genes, sample=SAMPLES)
    params: outdir = pjoin(OUT, "gene_tables"),
            samp = SAMPLES,
            samplelist = pjoin(OUT, "gene_tables", "amg_samplelist.txt"),
            workingdir = OUT,
            amg_heatmap = pjoin(OUT, "gene_tables", "amg_heatmap_cpm"),
            pythonpath = clust_conf["amg_tables"]["pythonpath"]
    output: amgs = pjoin(OUT, "gene_tables", "amg_cpm.tsv")
    shell:"""

    ## cleanup possible previous failed run
    rm -f {output.amgs}
    mkdir -p {params.outdir}

    ## make input sample list for script
    echo "{params.samp}" >{params.samplelist}
    tr " " "\\n" <{params.samplelist} >{params.outdir}/temp_amgs && mv {params.outdir}/temp_amgs {params.samplelist}


    ## dramv-distill amg_summary abund tables

echo "Collating abundances of AMGs." 

    python3 {config[scriptdir]}/scripts/dramv_amgs_table.py {params.workingdir} {params.samplelist} -v cpm >{output.amgs}

    ## heatmap of amgs
    export PYTHONPATH={params.pythonpath}
    python3 {config[scriptdir]}/scripts/plotnine_heatmap.py {output.amgs} {params.amg_heatmap} \
           -t "Heatmap of AMGs" -d "gene_description" -a "cpm" 

    rm {params.samplelist}

    """



rule iphop:
    threads: clust_conf["iphop"]["threads"]
    envmodules: *clust_conf["iphop"]["modules"]
    input: fasta = rules.mmseqs.output.renamed_DB_clu_rep_fasta
    params: outdir = pjoin(HOST_DIR)
    log: pjoin(HOST_DIR, "iphop.log")
    output: pjoin(HOST_DIR, "Host_prediction_to_genus_m75.csv")

    shell:"""
    ## iphop for bacteriophage host calls
    iphop --version


    ## cleanup possible previous failed run

    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    echo "Predicting viral hosts with iPHoP. This step may take a while. See log {log}." 

    iphop predict --fa_file {input.fasta} --db_dir {config[iphopdb]} --min_score {config[iphop_min_score]} \
	--out_dir {params.outdir} -t {threads} 1>>{log} 2>>{log}

    ## remove working dir
    rm -rf {params.outdir}/Wdir

    echo "iPHoP finished."

    """

rule iphop_abund:
    threads: clust_conf["iphop_abund"]["threads"]
    envmodules: *clust_conf["iphop_abund"]["modules"]
    input: votu_table = rules.votu.output.votu,
           host = rules.iphop.output
    params: outdir = pjoin(HOST_DIR, "tmp"),
            host_table = pjoin(HOST_DIR, "Host_cpm_table.tsv"),
    log: pjoin(HOST_DIR, "iphop.log")
    output: host_krona = pjoin(HOST_DIR, "Host.krona.html")
         
            

    shell:"""
  
    ## cleanup possible previous failed run

    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    # make host abundance table

    python3 {config[scriptdir]}/scripts/host_cpm.py {input.votu_table} {input.host} \
        > {params.host_table} 2>>{log}

    # make krona charts

    python3 {config[scriptdir]}/scripts/format_iphop_and_votu_for_krona.py {input.votu_table} \
        {input.host} {params.outdir} 1>>{log} 2>>{log}
       
    ktImportText -o {output.host_krona} {params.outdir}/*.txt 1>>{log} 2>>{log}
  

    """

## check if diamond db specified
if config["run_diamond"]:
    if ("diamonddb" not in config) or (config["diamonddb"] is None):
        raise ValueError(f"run_diamond value in config is 'yes', but diamonddb config value is not specified")
    else:
        DIAMOND_DB_NAME=os.path.splitext(os.path.basename(config["diamonddb"]))[0]
else:
     DIAMOND_DB_NAME=None

rule diamond:
    threads: clust_conf["diamond"]["threads"]
    envmodules: clust_conf["diamond"]["modules"]
    input: rules.genomad.output.proteins
    output: pjoin(SOUT, "diamond", "{sample}." + DIAMOND_DB_NAME + '.tsv') if config["run_diamond"] else "{sample}.temp.txt"
    log: pjoin(SOUT, "diamond", "{sample}" + ".diamond.log")
    params: outdir = pjoin(SOUT, "diamond"),
            s = "{sample}",
            tempdir = pjoin(clust_conf["diamond"]["tmpdir"], "{sample}_virome_diamond"),
            dbname = DIAMOND_DB_NAME
    shell:"""
    ## align gene sequences found by genomad to a diamond database for
    ## additional functional annotation
    diamond --version

    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir} 

    ## https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#memory--performance-options
    mkdir -p {params.tempdir}
    trap 'rm -rvf {params.tempdir}' EXIT

    echo "Annotating genes from geNomad for {wildcards.sample} against the nr database using diamond. See log file {log}." 

    diamond blastp --threads {threads} --max-target-seqs 2 -b 13 --tmpdir {params.tempdir} \
            --query {input} --db {config[diamonddb]} \
            --daa {params.outdir}/{params.s}.{params.dbname}.daa 1>>{log} 2>>{log}

    diamond view --threads {threads} --outfmt 6 qseqid pident qcovhsp scovhsp length mismatch gapopen qstart qend sstart send evalue bitscore stitle -a {params.outdir}/{params.s}.{params.dbname}.daa -o {output} 1>>{log} 2>>{log}

    ## add header to file
    sed -i '1s;^;qseqid\\tpident\\tqcovhsp\\tscovhsp\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tstitle\\n;' {output}

    echo "Diamond finished running for {wildcards.sample}."

    """



###### ALL RULE #############
rule all:
    input: GENOMADALL = expand(rules.genomad.output.fna, sample=SAMPLES),
           CHECKVALL = expand(rules.checkv_filter.output, sample=SAMPLES),
           VERSEGALL = expand(rules.abund_genomad.output.readcounts_virus, sample=SAMPLES),
           BBTOOLS_DEDUPEALL = rules.bbtools_dedupe.output.unique_seqs,
	   MMSEQSALL = rules.mmseqs.output.DB_clu_rep_fasta,
           VOTUALL = rules.votu.output.votu,
           IPHOPALL = rules.iphop.output if config["run_iphop"] else [],
           DIAMALL = expand(rules.diamond.output, sample=SAMPLES) if config["run_diamond"] else [],
           DRAMVALL = expand(rules.dramv.output, sample=SAMPLES),
           VS4DRAMVALL = expand(rules.vs4dramv.output, sample=SAMPLES) if config["run_amgs"] else [],
           AMGSALL = expand(rules.amgs.output, sample=SAMPLES) if config["run_amgs"] else [],
           VERSEDALL = expand(rules.abund_dramv.output.readcounts_genes, sample=SAMPLES),
           GENETABLESALL = rules.gene_tables.output.vogdb,
           VERSEAMGSALL = expand(rules.abund_amgs.output, sample=SAMPLES) if config["run_amgs"] else [],
           AMGTABLESALL = rules.amg_tables.output.amgs if config["run_amgs"] else [],
           IPHOPABUND = rules.iphop_abund.output
