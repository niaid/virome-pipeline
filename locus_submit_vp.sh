#!/usr/bin/env bash
#$ -N run_vp
#$ -cwd
#$ -j y
#$ -l h_vmem=8G
#$ -o .
#$ -m e
#$ -M your_address@niaid.nih.gov


module load snakemake/7.8.5-Python-3.10.5 || exit 1

## project config
configfile=/path/to/config.yaml
clusterfile=$(grep "clusterfile" ${configfile} | awk '{ print $2 }')

## output directory could avoid all this grep by just pasting the path into this sc
o1=$(grep "workdir" ${configfile} | awk '{ print $2 }')
o2=$(grep "outdir" ${configfile} | awk '{ print $2 }')
outputdir=$(realpath ${o1}/${o2})

## job submit command - can't use drmaa on locus :( - would have to ask them to install it...
clustercmd="qsub -l h_vmem={cluster.h_vmem} -j y -pe threaded {cluster.threads} {cluster.extra} -o ${outputdir}/{cluster.log}"

## remove/add --dryrun to test
snakemake -s /path/to/virome-pipeline/Snakefile --jobs 32 --jobname "{name}.{cluster.jobname}.{jobid}" \
	  --configfile ${configfile} \
	  --cluster "${clustercmd}" --cluster-config ${clusterfile} \
	  --nolock --keep-going --keep-incomplete --use-envmodules \
	  all
