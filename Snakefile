from os.path import join as pjoin
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
    output: pjoin(SOUT, "snakefake")
    params: outdir = SOUT
    shell:"""
    mkdir -p {params.outdir}
    touch {output}
    """



########## PIPELINE RULES ########

rule vs1:
    threads: clust_conf["vs1"]["threads"]
    envmodules: clust_conf["vs1"]["modules"]
    input: fake = ancient(rules.createsampledir.output),
           assembly = pjoin(IN, "{sample}" + config["assembly_suffix"])
    params: outdir = pjoin(SOUT, "vs1"),
            tempdir = pjoin(TMP, "vs1")
    output: pjoin(SOUT, 'vs1/final-viral-combined.fa')
    shell:"""
    ## virsorter first pass to id viral seqs in assemblies
    virsorter --version

    ## cleanup possible previous run
    rm -rf {params.outdir}

    ## do work in temp directory and move to out directory at the end
    mkdir -p {params.tempdir}
    trap 'mv -v {params.tempdir} {params.outdir}' ERR EXIT SIGTERM

    virsorter run --keep-original-seq -i {input.assembly} \
            -w {params.tempdir} --include-groups dsDNAphage,ssDNA \
            --min-length 1000 --min-score 0.5 -j {threads} all


    """

rule checkv:
    threads: clust_conf["checkv"]["threads"]
    envmodules: *clust_conf["checkv"]["modules"]
    input: rules.vs1.output
    params: outdir = pjoin(SOUT, "checkv"),
            tempdir = pjoin(TMP, "checkv")
    output: pjoin(SOUT, "checkv", "combined.fna")
    shell:"""
    ## run checkv to qc virsorter results and trim host regions left at the end of proviruses

    ## cleanup possible previous run
    rm -rf {params.outdir}

    ## do work in temp directory and move to out directory at the end
    mkdir -p {params.tempdir}
    trap 'mv -v {params.tempdir} {params.outdir}' ERR EXIT SIGTERM

    checkv end_to_end {input} {params.tempdir} -t {threads}

    cat {params.tempdir}/proviruses.fna {params.tempdir}/viruses.fna \
        >{params.tempdir}/combined.fna
    """

rule vs4dramv:
    threads: clust_conf["vs4dramv"]["threads"]
    envmodules: clust_conf["vs4dramv"]["modules"]
    input: rules.checkv.output
    params: outdir = pjoin(SOUT, "vs2"),
            tempdir = pjoin(TMP, "vs2")
    output: fasta = pjoin(SOUT, "vs2/for-dramv/final-viral-combined-for-dramv.fa"),
            tab = pjoin(SOUT, "vs2/for-dramv/viral-affi-contigs-for-dramv.tab")
    shell:"""
    ## run virsorter again to produce input for DRAM-v
    virsorter --version

    ## cleanup possible previous failed run
    rm -rf {params.outdir}

    ## do work in temp directory and move to out directory at the end
    mkdir -p {params.tempdir}
    trap 'mv -v {params.tempdir} {params.outdir}' ERR EXIT SIGTERM

    virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off \
        --prep-for-dramv -i {input} \
        -w {params.tempdir} --include-groups dsDNAphage,ssDNA --min-length 1000 \
        --min-score 0.5 -j {threads} all

    """


###### ALL RULE #############
rule all:
    input: VS1ALL = expand(rules.vs1.output, sample=SAMPLES),
           CHECKVALL = expand(rules.checkv.output, sample=SAMPLES),
           VS2ALL = expand(rules.vs4dramv.output.tab, sample=SAMPLES)
