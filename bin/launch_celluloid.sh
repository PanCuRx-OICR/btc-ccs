#!/usr/bin/env bash

# tumour id
tumor="$1"
normal_sample="$2"
donor="$3"
core_out_dir="$4"

# Check if is provided
if [ -z "$tumor" | -z "$normal_sample" | -z "$donor" ]; then
    echo "Usage: $0 <tumor> <normal_sample> <donor>"
    exit 1
fi

#connect resources
data_dir=/.mounts/labs/PCSI/references/GRCh38/dna/HMMcopy
script_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/celluloid/

#set input
input_from_local=/.mounts/labs/PCSI/pipeline/hg38_production
input_location=${input_from_local}/${donor}/${tumor}/wgs/bwa/0.7.17
normal_input_location=${input_from_local}/${donor}/${normal_sample}/wgs/bwa/0.7.17

#set output
out_dir=${core_out_dir}/${tumor}/celluloid/

### add two possible places for sex_file
#if [ -e "${input_location}/gendertype/final/${donor}.final_gender.txt" ]; then
#    sex_file=${input_location}/gendertype/final/${donor}.final_gender.txt
#else
#    sex_file=${input_location}/gendertype/${tumor}.xyr_gender.txt
#fi
sex_file=${core_out_dir}/${tumor}/sex/${normal_sample}.sex.txt


if [ -e "${out_dir}/${tumor}.celluloid.done" ]; then
      echo "file exists, skipping"
else

    mkdir ${out_dir}/${tumor}/
    touch ${out_dir}/${tumor}.celluloid.chk

    Rscript ${script_dir}/run_celluloidXY.R --tumorwig ${input_location}/HMMcopy/0.1.1/wigs/${tumor}.wig  --normalwig ${normal_input_location}/HMMcopy/0.1.1/wigs/${normal_sample}.wig --arfile ${input_location}/celluloidXY/v0.11.7/AR/${tumor}_AR.txt --output ${out_dir}/${tumor}  --gcwig ${data_dir}/GRCh38_gc.wig  --mapwig ${data_dir}/GRCh38_map.wig --sample.name ${tumor}  --sex_file ${sex_file} --script_dir ${script_dir} --colors new


    ln -s ${out_dir}/${tumor}/solution1 ${out_dir}/solution

    Rscript ${script_dir}/submit_celluloid_subclones.R --rda.dir ${out_dir}/${tumor}/saved_data --param.file ${out_dir}/solution/parameters_${tumor}.txt --sample.name ${tumor} --sex_file ${sex_file} --celluloid_dir ${script_dir} --output_path ${out_dir}/solution/

    gzip ${out_dir}/solution/segments_${tumor}.txt

    if [ -e "${out_dir}/solution/segments_${tumor}.txt.gz" ]; then
        rm   ${out_dir}/${tumor}.celluloid.chk
        touch ${out_dir}/${tumor}.celluloid.done
    fi
fi