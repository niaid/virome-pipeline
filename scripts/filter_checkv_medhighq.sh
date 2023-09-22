#!/usr/bin/env bash


# $1 is input checkv directory that has quality_summary
# $2 is output filename of filtered contigs that are only medium/high/complete (optional) default is ${1}/viruses_qualfilt.fna

defout=viruses_qualfilt.fna

if [ -z ${1+x} ]; then
    echo -e "Usage:"
    echo -e "$0 <input_checkv_dir> [optional output filename]"
    echo -e "if output file not given, will use input_checkv_dir/${defout}"
    exit 0
fi


contiglist=${1}/checkv.qualfilt.contigs.txt

if [ -z ${2+x} ]; then
    outfile=${1}/${defout}
else
    outfile=${2}
fi

module load seqtk

grep -e "Medium-quality" -e "High-quality" -e "Complete" $1/quality_summary.tsv | awk '{ print $1 }' >${contiglist}

seqtk subseq ${1}/combined.fna ${contiglist} > ${outfile}


rm -f ${contiglist}
