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

rule vs1:
    threads: clust_conf["vs1"]["threads"]
    envmodules: clust_conf["vs1"]["modules"]
    input: fake = ancient(rules.createsampledir.output),
           assembly = pjoin(IN, "{sample}" + config["assembly_suffix"])
    params: outdir = pjoin(SOUT, "vs1")
    output: pjoin(SOUT, 'vs1/final-viral-combined.fa')
    shell:"""
    ## virsorter first pass to id viral seqs in assemblies
    virsorter --version

    ## cleanup possible previous run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    virsorter run --keep-original-seq -i {input.assembly} \
            -w {params.outdir} --include-groups dsDNAphage,ssDNA \
            --min-length {config[vs_min_length]} --min-score 0.5 -j {threads} all


    """

rule checkv:
    threads: clust_conf["checkv"]["threads"]
    envmodules: *clust_conf["checkv"]["modules"]
    input: rules.vs1.output
    params: outdir = pjoin(SOUT, "checkv")
    output: pjoin(SOUT, "checkv", "combined.fna")
    shell:"""
    ## run checkv to qc virsorter results and trim host regions left at the end of proviruses

    ## cleanup possible previous run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    checkv end_to_end {input} {params.outdir} -d {config[checkvdb]} -t {threads}

    cat {params.outdir}/proviruses.fna {params.outdir}/viruses.fna \
        >{params.outdir}/combined.fna
    """

rule vs4dramv:
    threads: clust_conf["vs4dramv"]["threads"]
    envmodules: clust_conf["vs4dramv"]["modules"]
    input: rules.checkv.output
    params: outdir = pjoin(SOUT, "vs2")
    output: fasta = pjoin(SOUT, "vs2/for-dramv/final-viral-combined-for-dramv.fa"),
            tab = pjoin(SOUT, "vs2/for-dramv/viral-affi-contigs-for-dramv.tab")
    shell:"""
    ## run virsorter again to produce input for DRAM-v
    virsorter --version

    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off \
        --prep-for-dramv -i {input} \
        -w {params.outdir} --include-groups dsDNAphage,ssDNA --min-length {config[vs_min_length]} \
        --min-score 0.5 -j {threads} all

    """

rule dramv:
    threads: clust_conf["dramv"]["threads"]
    envmodules: *clust_conf["dramv"]["modules"]
    input: fasta = rules.vs4dramv.output.fasta,
           tab = rules.vs4dramv.output.tab
    params: outdir = pjoin(SOUT, "dramv")
    output: summary = pjoin(SOUT, "dramv/dramv-distill/amg_summary.tsv"),
            genes = pjoin(SOUT, "dramv/dramv-annotate/genes.faa")
    shell:"""
    ## DRAM-v annotation of viral sequences
    DRAM-setup.py version
    mmseqs -h | head -n 7

    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    DRAM-v.py annotate -i {input.fasta} \
        -v {input.tab} \
        -o {params.outdir}/dramv-annotate --skip_trnascan \
        --threads {threads} --min_contig_size {config[vs_min_length]}

    ## summarize annotations
    DRAM-v.py distill -i {params.outdir}/dramv-annotate/annotations.tsv \
       -o {params.outdir}/dramv-distill

    """


DIAMOND_DB_NAME=os.path.splitext(os.path.basename(config["diamonddb"]))[0]
rule diamond:
    threads: clust_conf["diamond"]["threads"]
    envmodules: clust_conf["diamond"]["modules"]
    input: rules.dramv.output.genes
    output: pjoin(SOUT, "diamond/{sample}." + DIAMOND_DB_NAME + '.tsv')
    params: outdir = pjoin(SOUT, "diamond"),
            s = "{sample}",
            tempdir = pjoin(clust_conf["diamond"]["tmpdir"], "{sample}_virome_diamond"),
            dbname = DIAMOND_DB_NAME
    shell:"""
    ## align gene sequences found by DRAM-v/prodigal to a diamond database for
    ## additional functional annotation
    diamond --version

    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    ## https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options#memory--performance-options
    mkdir -p {params.tempdir}
    trap 'rm -rvf {params.tempdir}' EXIT

    diamond blastp --threads {threads} --max-target-seqs 2 -b 13 --tmpdir {params.tempdir} \
            --query {input} --db {config[diamonddb]} \
            --daa {params.outdir}/{params.s}.{params.dbname}.daa

    diamond view --threads {threads} --outfmt 6 qseqid pident qcovhsp scovhsp length mismatch gapopen qstart qend sstart send evalue bitscore stitle -a {params.outdir}/{params.s}.{params.dbname}.daa -o {output}

    ## add header to file
    sed -i '1s;^;qseqid\\tpident\\tqcovhsp\\tscovhsp\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tstitle\\n;' {output}

    """



###### ALL RULE #############
rule all:
    input: VS1ALL = expand(rules.vs1.output, sample=SAMPLES),
           CHECKVALL = expand(rules.checkv.output, sample=SAMPLES),
           DRAMVALL = expand(rules.dramv.output.summary, sample=SAMPLES),
           DIAMALL = expand(rules.diamond.output, sample=SAMPLES)
