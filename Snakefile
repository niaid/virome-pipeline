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

rule genomad:
    threads: clust_conf["genomad"]["threads"]
    envmodules: clust_conf["genomad"]["modules"]
    input: fake = ancient(rules.createsampledir.output),
           assembly = pjoin(IN, "{sample}" + config["assembly_suffix"])
    params: outdir = pjoin(SOUT, "genomad"),
            renamed_assembly = pjoin(SOUT, "genomad", "{sample}" + ".fasta")
    output: fna = pjoin(SOUT, "genomad", "{sample}_summary", "{sample}_virus.fna"),
            gff = pjoin(SOUT, "genomad", "{sample}_summary", "{sample}_virus_genes.gff"),
            assembly_headermap = pjoin(SOUT, "genomad", "{sample}" + "_assembly_headermap.txt")
    shell:"""
    ## genomad search for viral contigs
    genomad --version

    ## cleanup possible previous run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

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
    genomad end-to-end --cleanup --threads {threads} $renamed_assembly {params.outdir} {config[genomaddb]}

    ## make gff from virus genes files
    python3 {config[scriptdir]}/scripts/genomad_genes2gff.py {params.outdir}/{wildcards.sample}_summary/{wildcards.sample}_virus_genes.tsv >{output.gff}

    rm $renamed_assembly

    """

rule verse:
    threads: clust_conf["verse"]["threads"]
    envmodules: *clust_conf["verse"]["modules"]
    input:  gff = rules.genomad.output.gff,
            bam = ancient(pjoin(IN, "{sample}" + config["bam_suffix"]))
    params: outdir = pjoin(SOUT, "verse"),
            prefix = pjoin(SOUT, "verse", "{sample}_virus_genes.count"),
            counts_only = pjoin(SOUT, "verse", "{sample}_virus_genes.count.CDS.txt")
    output: pjoin(SOUT, "verse", "{sample}_virus_genes.count.CDS.cpm.txt")
    shell:"""
    ## estimate gene abundances with verse
    verse -v

    ## cleanup possible previous run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    verse -a {input.gff} -o {params.prefix} -g ID -z 1 -t CDS -l -T {threads} {input.bam}

    python3 {config[scriptdir]}/scripts/calc_cpm.py {params.counts_only} >{output}

    """

rule checkv:
    threads: clust_conf["checkv"]["threads"]
    envmodules: *clust_conf["checkv"]["modules"]
    input: rules.genomad.output.fna
    params: outdir = pjoin(SOUT, "checkv")
    output: pjoin(SOUT, "checkv", "combined.fna")
    shell:"""
    ## run checkv to qc genomad results and trim host regions left at the end of proviruses

    ## cleanup possible previous run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    checkv end_to_end {input} {params.outdir} -d {config[checkvdb]} -t {threads}

    cat {params.outdir}/proviruses.fna {params.outdir}/viruses.fna \
        >{params.outdir}/combined.fna


    """

rule vsearch:
    threads: clust_conf["vsearch"]["threads"]
    envmodules: clust_conf["vsearch"]["modules"]
    input: expand(rules.genomad.output.fna, sample=SAMPLES)
    params: outdir = pjoin(OUT, "vsearch"),
            input_ctgs = pjoin(OUT, "vsearch", "all_input_contigs.fasta")
    output: info = pjoin(OUT, "vsearch", "info.txt"),
            consensus = pjoin(OUT, "vsearch", "consensus_sequences.fasta")
    shell:"""
    ## run vsearch on all contigs to remove duplicates

    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    cat {input} >{params.input_ctgs}

    vsearch --cluster_size {params.input_ctgs} --id 0.95 --consout {output.consensus} \
      --clusterout_id --maxseqlength 500000 --threads 16 --iddef 0 --minseqlength 5000 \
      --uc {output.info}


    ## cleanup remove contigs to save space
    rm {params.input_ctgs}

    """


rule dramv:
    threads: clust_conf["dramv"]["threads"]
    envmodules: *clust_conf["dramv"]["modules"]
    input: fasta = rules.genomad.output.fna
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
    input: fasta = rules.genomad.output.fna
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
    input: GENOMADALL = expand(rules.genomad.output.gff, sample=SAMPLES),
           CHECKVALL = expand(rules.checkv.output, sample=SAMPLES),
           DRAMVALL = expand(rules.dramv.output, sample=SAMPLES),
	   IPHOPALL = expand(rules.iphop.output, sample=SAMPLES),
           DIAMALL = expand(rules.diamond.output, sample=SAMPLES),
           VERSEALL = expand(rules.verse.output, sample=SAMPLES),
           VSEARCHALL = rules.vsearch.output.info
