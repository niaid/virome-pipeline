from os.path import join as pjoin
import os
import snakemake.utils as sm

### cluster config
scriptdir = config["scriptdir"] ## virome-pipeline repo/script directory
include: pjoin(scriptdir, "scripts", "cluster_setup.smk")

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

rule genomad:
    threads: clust_conf["genomad"]["threads"]
    envmodules: clust_conf["genomad"]["modules"]
    input: fake = ancient(rules.createsampledir.output),
           assembly = pjoin(IN, "{sample}" + config["assembly_suffix"])
    params: outdir = pjoin(SOUT, "genomad"),
            genomad2gff = pjoin(scriptdir, "scripts", "genomad_genes2gff.py"),
            python_exe = clust_conf["genomad"]["python_exe"]
    output: pjoin(SOUT, "genomad", "{sample}_summary", "{sample}_virus.fna")
    shell:"""
    ## genomad search for viral contigs
    genomad --version

    ## cleanup possible previous run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    ## link assembly to name that works better for genomad
    ln -s {input.assembly} {params.outdir}/{sample}.fasta

    genomad end-to-end --cleanup {params.outdir}/{sample}.fasta {params.outdir} {config[genomaddb]}

    ## make gff from final genes files
    {params.python_exe} {params.genomad2gff} {params.outdir}/{sample}_summary/{sample}_virus_genes.tsv >{params.outdir}/{sample}_summary/{sample}_virus_genes.gff

    rm {params.outdir}/{sample}.fasta

    """

rule checkv:
    threads: clust_conf["checkv"]["threads"]
    envmodules: *clust_conf["checkv"]["modules"]
    input: rules.genomad.output
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

rule dramv:
    threads: clust_conf["dramv"]["threads"]
    envmodules: *clust_conf["dramv"]["modules"]
    input: fasta = rules.genomad.output
    params: outdir = pjoin(SOUT, "dramv")
    output: genes = pjoin(SOUT, "dramv/dramv-annotate/genes.faa")

    shell:"""
    ## DRAM-v annotation of viral sequences
    DRAM-setup.py version
    mmseqs -h | head -n 7

    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    DRAM-v.py annotate -i {input.fasta} \
        -o {params.outdir}/dramv-annotate --skip_trnascan \
        --threads {threads} --min_contig_size {config[vs_min_length]}
    """

rule iphop:
    threads: clust_conf["iphop"]["threads"]
    envmodules: *clust_conf["iphop"]["modules"]
    input: fasta = rules.genomad.output
    params: outdir = pjoin(SOUT, "iphop")
    output: genes = pjoin(SOUT, "iphop/Host_prediction_to_genome_m90.csv")

    shell:"""
    ## iphop for bacteriophage host calls
    iphop --version


    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}


iphop predict --fa_file {input.fasta} --db_dir {config[iphopdb]} \
	--out_dir {params.outdir} -t {threads}


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
    input: GENOMADALL = expand(rules.genomad.output, sample=SAMPLES),
           CHECKVALL = expand(rules.checkv.output, sample=SAMPLES),
           DRAMVALL = expand(rules.dramv.output, sample=SAMPLES),
	   IPHOPALL = expand(rules.iphop.output, sample=SAMPLES),
           DIAMALL = expand(rules.diamond.output, sample=SAMPLES)
