#!/usr/bin/env bash
#$ -N run_vp
#$ -cwd
#$ -j y
#$ -l h_vmem=8G
#$ -o .
#$ -m e
#$ -M your_address@niaid.nih.gov


module load snakemake || exit 1
## project config
configfile=/path/to/config.yaml
clusterfile=$(grep "clusterfile" ${configfile} | awk '{ print $2 }')

## output directory could avoid all this grep by just hard coding it
o1=$(grep "workdir" ${configfile} | awk '{ print $2 }')
o2=$(grep "outdir" ${configfile} | awk '{ print $2 }')
outputdir=$(realpath ${o1}/${o2})

## job submit command
drmaacmd=" -l h_vmem={cluster.h_vmem} -j y -pe threaded {cluster.threads} {cluster.extra}"


snakemake -s /path/to/virome-pipeline/Snakefile --jobs 32 --jobname "{name}.{cluster.jobname}.{jobid}" \
	  --configfile ${configfile} \
	  --drmaa "${drmaacmd}" --cluster-config ${clusterfile} --drmaa-log-dir ${outputdir}/{cluster.log} \
	  --nolock --keep-going --keep-incomplete --use-envmodules
