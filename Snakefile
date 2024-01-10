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
    output: fna = pjoin(SOUT, "genomad", "{sample}_summary", "{sample}_virus.fna"),
            genes = pjoin(SOUT, "genomad", "{sample}_summary", "{sample}_virus_genes.tsv"),
            proteins =  pjoin(SOUT, "genomad", "{sample}_summary", "{sample}_virus_proteins.faa"),
            summary = pjoin(SOUT, "genomad", "{sample}_summary","{sample}_virus_summary.tsv"),
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
    genomad end-to-end {params.genomad_filter} --enable-score-calibration --cleanup --threads {threads} $renamed_assembly {params.outdir} {config[genomaddb]}

    rm $renamed_assembly

    """

rule verse_genomad:
    threads: clust_conf["verse"]["threads"]
    envmodules: *clust_conf["verse"]["modules"]
    input:  genes = rules.genomad.output.genes,
            fna = rules.genomad.output.fna,
            virus = rules.genomad.output.summary,
            headermap = rules.genomad.output.assembly_headermap,
            bam = ancient(pjoin(IN, "{sample}" + config["bam_suffix"]))
    params: outdir = pjoin(SOUT, "verse_genomad"),
            prefix_genes = pjoin(SOUT, "verse_genomad", "{sample}_virus_genes.count"),
            counts_only_genes = pjoin(SOUT, "verse_genomad", "{sample}_virus_genes.count.CDS.txt"),
            prefix_virus = pjoin(SOUT, "verse_genomad", "{sample}_virus.count"),
            counts_only_virus = pjoin(SOUT, "verse_genomad", "{sample}_virus.count.CDS.txt")
    output: readcounts_genes = pjoin(SOUT, "verse_genomad", "{sample}_virus_genes.count.CDS.cpm.txt"),
            readcounts_virus = pjoin(SOUT, "verse_genomad", "{sample}_virus.count.CDS.cpm.txt")
    shell:"""
    ## estimate gene abundances with verse
    verse -v

    ## cleanup possible previous run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}


    ## make gff from virus genes files
    python3 {config[scriptdir]}/scripts/genomad_genes2gff.py {input.genes} -m {input.headermap} >{params.prefix_genes}.gff
    verse -a {params.prefix_genes}.gff -o {params.prefix_genes} -g ID -z 1 -t CDS -l -T {threads} {input.bam}

    python3 {config[scriptdir]}/scripts/calc_cpm.py {params.counts_only_genes} >{output.readcounts_genes}

    ## make gff from virus summary; use default -z 1 instead of -z 5
    python3 {config[scriptdir]}/scripts/genomad_virus2gff.py {input.virus} -m {input.headermap} >{params.prefix_virus}.gff
    verse -a {params.prefix_virus}.gff -o {params.prefix_virus} -g ID -z 1 -t CDS -l -T {threads} {input.bam}
    python3 {config[scriptdir]}/scripts/calc_cpm.py {params.counts_only_virus} >{output.readcounts_virus}

    """


rule checkv:
    threads: clust_conf["checkv"]["threads"]
    envmodules: *clust_conf["checkv"]["modules"]
    input: rules.genomad.output.fna
    params: outdir = pjoin(SOUT, "checkv")
    output: fna = pjoin(SOUT, "checkv", "combined.fna"),
            qsum = pjoin(SOUT, "checkv", "quality_summary.tsv")
    shell:"""
    ## run checkv to qc genomad results and trim host regions left at the end of proviruses

    ## cleanup possible previous run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    checkv end_to_end {input} {params.outdir} -d {config[checkvdb]} -t {threads}

    cat {params.outdir}/proviruses.fna {params.outdir}/viruses.fna >{output.fna}

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
            tempclist = pjoin(SOUT, "checkv", "tempclist.txt")
    output: pjoin(SOUT, "checkv", "checkv_filtered_genomad_viruses.fna")
    shell:"""
    ## use checkv quality assessment to filter viral contigs for clustering

    ## cleanup possible previous run
    rm -rf {output}


    if [ "{params.checkv_q}" = "complete" ]; then
       grep -e "Complete" {input.qsum} | awk '{{ print $1 }}' >{params.tempclist}
    elif [ "{params.checkv_q}" = "high" ]; then
       grep -e "Complete" -e "High-quality" {input.qsum} | awk '{{ print $1 }}' >{params.tempclist}
    elif [ "{params.checkv_q}" = "medium" ]; then
       grep -e "Complete" -e "High-quality" -e "Medium-quality" {input.qsum} | awk '{{ print $1 }}' >{params.tempclist}
    fi

    if [ "{params.checkv_q}" != "all" ]; then
        seqtk subseq {input.fna} {params.tempclist} > {params.outdir}/temp && mv {params.outdir}/temp {output}
        rm -fv {params.tempclist}
    else
        cp {input.fna} {output}
    fi

    """



rule bbmap:
    threads: clust_conf["bbmap"]["threads"]
    envmodules: clust_conf["bbmap"]["modules"]
    input: expand(rules.checkv_filter.output, sample=SAMPLES)
    params: outdir = pjoin(OUT, "bbmap"),
            input_ctgs = pjoin(OUT, "bbmap", "all_input_contigs.fasta"),
	    log = pjoin(OUT, "bbmap", "log.txt"),
	    stats = pjoin(OUT, "bbmap", "cluster_stats.txt")
    output: unique_seqs = pjoin(OUT, "bbmap", "unique_seqs.fasta")
    shell:"""
    ## run bbmap on all contigs to remove duplicates
    dedupe.sh --version

    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    cat {input} >{params.input_ctgs}

    dedupe.sh in={params.input_ctgs} out={output.unique_seqs} csf={params.stats} minscaf={config[vs_min_length]} \
	mergenames=t ex=f usejni=t threads={threads}


    """

rule mmseqs:
    threads: clust_conf["mmseqs"]["threads"]
    envmodules: *clust_conf["mmseqs"]["modules"]
    input: rules.bbmap.output.unique_seqs
    params: outdir = pjoin(OUT, "mmseqs"),
            DB_dir = pjoin(OUT, "mmseqs", "DB"),
            DB = pjoin(OUT, "mmseqs", "DB/DB"),
	    DB_clu = pjoin(OUT, "mmseqs", "DB/DB_clu"),
	    DB_clu_seq = pjoin(OUT, "mmseqs", "DB/DB_clu_seq"),
	    DB_clu_fasta = pjoin(OUT, "mmseqs", "cluster_seqs.fasta"),
	    DB_clu_rep = pjoin(OUT, "mmseqs", "DB/DB_clu_rep")
    output: DB_clu_rep_fasta = pjoin(OUT, "mmseqs", "representative_seqs.fasta"),
            DB_clu_tsv = pjoin(OUT, "mmseqs", "DB_clu.tsv"),
            flat_DB_clu_tsv = pjoin(OUT, "mmseqs", "flat_DB_clu.tsv"),
            renamed_DB_clu_rep_fasta = pjoin(OUT, "mmseqs", "representative_seqs.renamed.fasta")
    shell:"""
    ## run mmseqs on all deduped genomes
    mmseqs version

    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    mkdir {params.DB_dir}

    mmseqs createdb {input} {params.DB}

    mmseqs cluster {params.DB} {params.DB_clu} $TMPDIR --cov-mode 1 -c 0.85 \
	--min-seq-id 0.95 --cluster-mode 2 --threads {threads}

    mmseqs createtsv {params.DB} {params.DB} {params.DB_clu} {output.DB_clu_tsv} --threads {threads}

    mmseqs createseqfiledb {params.DB} {params.DB_clu} {params.DB_clu_seq} --threads {threads}

    mmseqs result2flat {params.DB} {params.DB} {params.DB_clu_seq} {params.DB_clu_fasta}

    mmseqs createsubdb {params.DB_clu} {params.DB} {params.DB_clu_rep}

    mmseqs convert2fasta {params.DB_clu_rep} {output.DB_clu_rep_fasta}

    ## rename repseq to only use first dup sequence from bbmap
    # flatten clu.tsv file so each row is one repseq and one seq
    python3 {config[scriptdir]}/scripts/flatten_mmseqs_tsv.py {output.DB_clu_tsv} | sort -u > {output.flat_DB_clu_tsv}
    # rename representative_seqs.fasta
    sed '/^>/ s!>[^>]*!!2g' {output.DB_clu_rep_fasta} >{output.renamed_DB_clu_rep_fasta}


    """

rule votu:
    threads: clust_conf["votu"]["threads"]
    envmodules: clust_conf["votu"]["modules"]
    input: mmseqs = rules.mmseqs.output.flat_DB_clu_tsv,
           abund = expand(rules.verse_genomad.output.readcounts_virus, sample=SAMPLES),
           summ = expand(rules.genomad.output.summary, sample=SAMPLES)
    params: outdir = pjoin(OUT, "votu"),
            filelist = pjoin(OUT, "votu", "abundfiles.txt"),
            summlist = pjoin(OUT, "votu", "summaryfiles.txt")
    output: votu = pjoin(OUT, "votu", "vOTU_table.tsv"),
            gmdanno = pjoin(OUT, "votu", "repseq_genomad_virus_summary.tsv")
    shell:"""
    ## make vOTU table
    python3 --version

    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    ## vOTU table
    echo "{input.abund}" >{params.filelist}
    tr " " "\\n" <{params.filelist} >{params.outdir}/temp && mv {params.outdir}/temp {params.filelist}
    python3 {config[scriptdir]}/scripts/make_votu_table.py {input.mmseqs} {params.filelist} >{output.votu}

    ## collate genomad annotations
    echo "{input.summ}" | tr " " "\\n" >{params.summlist}
    python3 {config[scriptdir]}/scripts/repseq_genomad_annotations.py {output.votu} {params.summlist} >{output.gmdanno}

    """

rule vs4dramv:
    threads: clust_conf["vs4dramv"]["threads"]
    envmodules: clust_conf["vs4dramv"]["modules"]
    input: rules.checkv.output.fna
    params: outdir = pjoin(SOUT, "vs2")
    output: fasta = pjoin(SOUT, "vs2/for-dramv/final-viral-combined-for-dramv.fa"),
            tab = pjoin(SOUT, "vs2/for-dramv/viral-affi-contigs-for-dramv.tab")
    shell:"""
    ## run virsorter to produce input for DRAM-v
    virsorter --version

    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    virsorter run --seqname-suffix-off --viral-gene-enrich-off --provirus-off \
        --prep-for-dramv -i {input} -w {params.outdir} --include-groups dsDNAphage,ssDNA,NCLDV \
        --min-length {config[vs_min_length]} --min-score 0.5 -j {threads} all

    """

rule amgs:
    threads: clust_conf["amgs"]["threads"]
    envmodules: clust_conf["amgs"]["modules"]
    input: fasta = rules.vs4dramv.output.fasta,
           tab = pjoin(SOUT, "vs2/for-dramv/viral-affi-contigs-for-dramv.tab")
    params: outdir = pjoin(SOUT, "amgs")
    output: genes = pjoin(SOUT, "amgs/dramv-annotate/genes.faa")

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
       
        DRAM-v.py distill -i {params.outdir}/dramv-annotate/annotations.tsv \
        -o {params.outdir}/dramv-distill


    """

rule dramv:
    threads: clust_conf["dramv"]["threads"]
    envmodules: clust_conf["dramv"]["modules"]
    input: fasta = rules.checkv.output.fna
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
    input: fasta = rules.mmseqs.output.renamed_DB_clu_rep_fasta
    params: outdir = pjoin(OUT, "iphop")
    output: genes = pjoin(OUT, "iphop", "Host_prediction_to_genome_m90.csv")

    shell:"""
    ## iphop for bacteriophage host calls
    iphop --version


    ## cleanup possible previous failed run
    rm -rf {params.outdir}
    mkdir -p {params.outdir}

    iphop predict --fa_file {input.fasta} --db_dir {config[iphopdb]} \
	--out_dir {params.outdir} -t {threads}

    """

if "diamonddb" in config:
    DIAMOND_DB_NAME=os.path.splitext(os.path.basename(config["diamonddb"]))[0]
else:
    DIAMOND_DB_NAME=None

rule diamond:
    threads: clust_conf["diamond"]["threads"]
    envmodules: clust_conf["diamond"]["modules"]
    input: rules.genomad.output.proteins
    output: pjoin(SOUT, "diamond", "{sample}." + DIAMOND_DB_NAME + '.tsv') if DIAMOND_DB_NAME else "temp.txt"
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

    diamond blastp --threads {threads} --max-target-seqs 2 -b 13 --tmpdir {params.tempdir} \
            --query {input} --db {config[diamonddb]} \
            --daa {params.outdir}/{params.s}.{params.dbname}.daa

    diamond view --threads {threads} --outfmt 6 qseqid pident qcovhsp scovhsp length mismatch gapopen qstart qend sstart send evalue bitscore stitle -a {params.outdir}/{params.s}.{params.dbname}.daa -o {output}

    ## add header to file
    sed -i '1s;^;qseqid\\tpident\\tqcovhsp\\tscovhsp\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tstitle\\n;' {output}

    """



###### ALL RULE #############
rule all:
    input: GENOMADALL = expand(rules.genomad.output.fna, sample=SAMPLES),
           CHECKVALL = expand(rules.checkv_filter.output, sample=SAMPLES),
           VERSEGALL = expand(rules.verse_genomad.output.readcounts_virus, sample=SAMPLES),
           BBMAPALL = rules.bbmap.output.unique_seqs,
	   MMSEQSALL = rules.mmseqs.output.DB_clu_rep_fasta,
           VOTUALL = rules.votu.output.votu,
           IPHOPALL = rules.iphop.output,
           DIAMALL = expand(rules.diamond.output, sample=SAMPLES) if DIAMOND_DB_NAME else [],
           DRAMVALL = expand(rules.dramv.output, sample=SAMPLES),
           VS4DRAMVALL = expand(rules.vs4dramv.output, sample=SAMPLES),
           AMGSALL = expand(rules.amgs.output, sample=SAMPLES)