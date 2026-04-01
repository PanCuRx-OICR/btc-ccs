#!/usr/bin/env bash

sample="$1"
donor="$2"
dir="$3"
base_dir="$4"

# Check if is provided
if [ -z "$sample" ]; then
    echo "Usage: $0 <sample> <donor> <directory>"
    exit 1
fi

module load samtools/1.17

sex_dir=${base_dir}/${dir}/sex/
cd $sex_dir

this_cram=/.mounts/labs/PCSI/pipeline/hg38_production/${donor}/${sample}/wgs/bwa/0.7.17/collapsed/${sample}.cram
this_bam=/.mounts/labs/PCSI/pipeline/hg38_production/${donor}/${sample}/wgs/bwa/0.7.17/collapsed/${sample}.bam

if [ -e "${sex_dir}/${sample}.sex.done" ]; then
    echo "file exists, skipping"
else
    touch ${sex_dir}/${sample}.sex.chk

    if test -f $this_cram
    then
        echo "viewing ${this_cram}"
        samtools view $this_cram  -F 2820 chrX chrY | cut -f3 | uniq -c | awk -v sample=$sample '{print sample"\t"$1"\t"$2}' >${sex_dir}/${sample}.XY.coverage.txt
    else
        echo "viewing ${this_bam}"
        samtools view $this_bam  -F 2820 chrX chrY | cut -f3 | uniq -c | awk -v sample=$sample '{print sample"\t"$1"\t"$2}' >${sex_dir}/${sample}.XY.coverage.txt
    fi

    echo "calculating ratios"

    xreads=$(awk 'NR==1 {print $2}' ${sex_dir}/${sample}.XY.coverage.txt)
    yreads=$(awk 'NR==2 {print $2}' ${sex_dir}/${sample}.XY.coverage.txt)

    # Divide the numbers
    y_to_x_reads=$(echo "$yreads / $xreads" | bc -l)

    # Take the natural logarithm of the result
    log_result=$(echo "-l($y_to_x_reads)" | bc -l)

    if (( $(echo "$log_result < 2.5" | bc -l) )); then
        echo "XY" > ${sex_dir}/${sample}.sex.txt
    fi 

    if (( $(echo "$log_result > 2.5" | bc -l) )); then
        echo "XX" > ${sex_dir}/${sample}.sex.txt
    fi

    if [ -e "${sex_dir}/${sample}.sex.txt" ]; then
        rm   ${sex_dir}/${sample}.sex.chk
        touch ${sex_dir}/${sample}.sex.done
    fi
fi