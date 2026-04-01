#!/usr/bin/env bash

# tumour id
tumour="$1"
donor="$2"


# Check if is provided
if [ -z "$tumour" ]; then
    echo "Usage: $0 <tumour> <donor>"
    exit 1
fi

work_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/${tumour}/amplicon/
cd $work_dir

input_location=/.mounts/labs/PCSI/pipeline/hg38_production/${donor}/${tumour}/wgs/bwa/0.7.17

tumor_cram=${input_location}/collapsed/${tumour}.cram
tumor_bam=${input_location}/collapsed/${tumour}.bam


if [ -e "${work_dir}/${tumour}.amplicon.done" ]; then
      echo "file exists, skipping"
else

    touch ${work_dir}/${tumour}.amplicon.chk

    mkdir ${work_dir}/results/
    cd ${work_dir}/results/
    rm -r ${work_dir}/results/*

    if test -f $tumor_cram
    then

    AmpliconSuite-pipeline.py -s ${tumour} -t 12  --bam $tumor_cram --run_AA --run_AC --ref GRCh38 --cngain 8 --downsample 8

    else
    AmpliconSuite-pipeline.py -s ${tumour} -t 12  --bam $tumor_bam --run_AA --run_AC --ref GRCh38 --cngain 8 --downsample 8 
    
    fi

    if [ -e "${work_dir}/results/${tumour}_classification/${tumour}_gene_list.tsv" ]; then
            rm   ${work_dir}/${tumour}.amplicon.chk
            touch ${work_dir}/${tumour}.amplicon.done
    fi
fi
