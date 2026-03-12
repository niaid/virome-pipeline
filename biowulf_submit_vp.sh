#!/usr/bin/env bash
### Biowulf params
#SBATCH --job-name=run_sm
#SBATCH --output=%x-%j.out
#SBATCH --partition=norm
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --cpus-per-task=1              # CPUs
#SBATCH --mem=8G                        # Job memory request
#SBATCH --time=2-00:00:00
#SBATCH --export=NONE
#SBATCH --gres=lscratch:1



module load snakemake/7.32.4 || exit 1

## protect us from user libs
export PYTHONNOUSERSITE=True

## project config
configfile=/data/subramanianp4/test/project_config.yaml

## cluster config file
clusterfile=$(grep "clusterfile" ${configfile} | awk '{ print $2 }')

## output directory could avoid all this grep by just pasting the path into this script
o1=$(grep "workdir" ${configfile} | awk '{ print $2 }')
o2=$(grep "outdir" ${configfile} | awk '{ print $2 }')
outputdir=$(realpath ${o1}/${o2})
mkdir -p ${outputdir}

## script dir
scriptdir=$(grep "scriptdir" ${configfile} | awk '{ print $2 }')

## job submit command - can't use drmaa on locus :( - would have to ask them to install it...
sbatchcmd="sbatch -c {cluster.threads} --mem={cluster.mem} --output=${outputdir}/{cluster.log}/%x-%j.out --partition={cluster.partition} --time={cluster.time} {cluster.extra}"


## add --dryrun to test
## sometimes useful to add --ignore-incomplete if snakemake's incomplete checking is not working for you.
## and --rerun-triggers mtime if you are debugging and don't want to keep re-running rules when code/format is changed but output is unchanged
snakemake -s ${scriptdir}/Snakefile --jobs 32 --jobname "{name}.{cluster.jobname}.{jobid}" \
	  --configfile ${configfile} \
	  --cluster "${sbatchcmd}" \
	  --cluster-config ${clusterfile} \
	  --cluster-cancel "scancel" \
	  --latency-wait 120 --max-jobs-per-second 1 \
	  --nolock --keep-going --keep-incomplete --use-envmodules \
	  all
