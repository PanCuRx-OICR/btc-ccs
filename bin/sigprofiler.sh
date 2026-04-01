#!/usr/bin/env bash

# tumour id
#put in argeuments instead
tumour_id="$1"
donor="$2"

# Check if is provided
if [ -z "$tumour_id" ]; then
    echo "Usage: $0 <tumour_id>"
    exit 1
fi

script_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/btc.scripts/

input_location=/.mounts/labs/PCSI/pipeline/hg38_production/${donor}/${tumour_id}/wgs/bwa/0.7.17
celluloid_location=/.mounts/labs/PCSI/users/fbeaudry/more.data/${tumour_id}/celluloid/solution/

outpath=/.mounts/labs/PCSI/users/fbeaudry/more.data/${tumour_id}/sigprofiler

SBS_signatures_path=${script_dir}/COSMIC_SBS96_Signatures.txt
DINUC_signatures_path=${script_dir}/COSMIC_DBS78_Signatures.txt
INDEL_signatures_path=${script_dir}/COSMIC_ID83_Signatures.txt

if [ -e "${outpath}/96/${tumour_id}.sigprofiler.done" ]; then
      echo "file exists, skipping"
else
    
    rm -r ${outpath}/96/*

    mkdir -p ${outpath}/96/input/

    touch ${outpath}/96/${tumour_id}.sigprofiler.chk

    python3 ${script_dir}/parse_SNV_freqs.py ${input_location}/final_strelka2_mutect2/${tumour_id}_merged_final.vcf.gz ${outpath}/96/${tumour_id} ${celluloid_location}/parameters_${tumour_id}.txt ${celluloid_location}/segments_${tumour_id}.txt.gz

    barebones_path=${outpath}/96/${tumour_id}.snv.barebones.vcf
    sed 's/chr//g'  ${barebones_path}  | awk -v tumour_id=$tumour_id '$1 !~ "pos1" {print $1"\t"$2"\t"tumour_id"\t"$4"\t"$5}'  >${outpath}/96/input/${tumour_id}.snv.cosmic.vcf

    #sigprofiler will try to read anything in the input directory, so make sure you've only got vcf files in there
    python3 ${script_dir}/sigProfilerAssignment.py ${outpath}/96/input/ ${outpath}/96/ ${tumour_id} ${SBS_signatures_path} 96

    if [ -e "${outpath}/96/${tumour_id}.most_probable_signatures.txt" ]; then
            rm   ${outpath}/96/${tumour_id}.sigprofiler.chk
            touch ${outpath}/96/${tumour_id}.sigprofiler.done
    fi
fi

if [ -e "${outpath}/DINUC/${tumour_id}.sigprofiler.done" ]; then
      echo "file exists, skipping"
else
    
    rm -r ${outpath}/DINUC/*

    mkdir -p ${outpath}/DINUC/input/

    touch ${outpath}/DINUC/${tumour_id}.sigprofiler.chk


    barebones_path=${outpath}/96/${tumour_id}.dbs.barebones.vcf
    sed 's/chr//g'  ${barebones_path}  | awk -v tumour_id=$tumour_id '$1 !~ "pos1" {print $1"\t"$2"\t"tumour_id"\t"$4"\t"$5}'  >${outpath}/DINUC/input/${tumour_id}.dbs.cosmic.vcf

    #sigprofiler will try to read anything in the input directory, so make sure you've only got vcf files in there
    python3 ${script_dir}/sigProfilerAssignment.py ${outpath}/DINUC/input/ ${outpath}/DINUC/ ${tumour_id} ${DINUC_signatures_path} DINUC

    if [ -e "${outpath}/DINUC/${tumour_id}.most_probable_signatures.txt" ]; then
            rm   ${outpath}/DINUC/${tumour_id}.sigprofiler.chk
            touch ${outpath}/DINUC/${tumour_id}.sigprofiler.done
    fi
fi

if [ -e "${outpath}/ID/${tumour_id}.sigprofiler.done" ]; then
      echo "file exists, skipping"
else
    
    rm -r ${outpath}/ID/*

    mkdir -p ${outpath}/ID/input/

    touch ${outpath}/ID/${tumour_id}.sigprofiler.chk

    if test -f ${input_location}/final_indel/${tumour_id}_merged_final_indels.vcf.gz
    then

    #output of this script is a directory + the root file name
    python3  ${script_dir}/parse_indel_freqs.py ${input_location}/final_indel/${tumour_id}_merged_final_indels.vcf.gz ${outpath}/ID/${tumour_id} ${celluloid_location}/parameters_${tumour_id}.txt ${celluloid_location}/segments_${tumour_id}.txt.gz

    else

    python3  ${script_dir}/parse_indel_freqs.py ${input_location}/final_indel.old/${tumour_id}_merged_final_indels.vcf.gz ${outpath}/ID/${tumour_id} ${celluloid_location}/parameters_${tumour_id}.txt ${celluloid_location}/segments_${tumour_id}.txt.gz

    fi

    barebones_path=${outpath}/ID/${tumour_id}.indel.barebones.vcf
    sed 's/chr//g'  ${barebones_path}  | awk -v tumour_id=$tumour_id '$1 !~ "pos1" {print $1"\t"$2"\t"tumour_id"\t"$4"\t"$5}'  >${outpath}/ID/input/${tumour_id}.indel.cosmic.vcf

    #sigprofiler will try to read anything in the input directory, so make sure you've only got vcf files in there
    python3 ${script_dir}/sigProfilerAssignment.py ${outpath}/ID/input/ ${outpath}/ID/ ${tumour_id} ${INDEL_signatures_path} ID

    if [ -e "${outpath}/ID/${tumour_id}.most_probable_signatures.txt" ]; then
            rm   ${outpath}/ID/${tumour_id}.sigprofiler.chk
            touch ${outpath}/ID/${tumour_id}.sigprofiler.done
    fi
fi

rm /.mounts/labs/PCSI/users/fbeaudry/more.data/${tumour_id}/sigprofiler/${tumour_id}.all.most_probable_signatures.txt
for this_type in 96 DINUC ID
do
cat /.mounts/labs/PCSI/users/fbeaudry/more.data/${tumour_id}/sigprofiler/${this_type}/${tumour_id}.most_probable_signatures.txt >>/.mounts/labs/PCSI/users/fbeaudry/more.data/${tumour_id}/sigprofiler/${tumour_id}.all.most_probable_signatures.txt
done
