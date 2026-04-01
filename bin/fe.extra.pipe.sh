#!/usr/bin/env bash

# tumour id
sample="$1"
normal="$2"
donor="$3"

# Check if is provided
if [ -z "$sample" ]; then
    echo "Usage: $0 <sample> <normal> <donor>"
    exit 1
fi

#sample=BTC_0022_Lv_M_526
#normal=BTC_0022_Ly_R
#donor=BTC0022

script_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/btc.scripts/
work_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/${sample}/
cd  ${work_dir}


#### run SV analyses
#conda install -c conda-forge r-base=4.2.0
#conda install -c conda-forge r-gmp r-optparse r-tidyr r-dplyr r-jsonlite r-plyr r-data.table r-mass r-gplots  r-fields  r-gensa r-cowplot r-devtools
#conda install bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg38
#conda install bioconda::bioconductor-structuralvariantannotation
#conda install bioconda::bioconductor-hmmcopy

#### sv

if [ -e "${work_dir}/sv.strict/${sample}.sv_strict.done" ]; then
      echo "file exists, skipping"
else
    
    mkdir ${work_dir}/sv.strict/
    touch ${work_dir}/sv.strict/${sample}.sv_strict.chk
    Rscript ${script_dir}/join_sv.R -t ${sample} -d ${donor} --output ${work_dir}/sv.strict/

    if [ -e "${work_dir}/sv.strict/${sample}.sv.cnv.txt" ]; then
            rm   ${work_dir}/sv.strict/${sample}.sv_strict.chk
            touch ${work_dir}/sv.strict/${sample}.sv_strict.done
    fi
fi

##hrdetect ##

input_location=/.mounts/labs/PCSI/pipeline/hg38_production/${donor}/${sample}/wgs/bwa/0.7.17

indel_vcf_location=${input_location}/final_indel/${sample}_merged_final_indels.vcf.gz

if test -f ${input_location}/final_indel/${sample}_merged_final_indels.vcf.gz
then
indel_vcf_location=${input_location}/final_indel/${sample}_merged_final_indels.vcf.gz
else
indel_vcf_location=${input_location}/final_indel.old/${sample}_merged_final_indels.vcf.gz
fi

LOH_seg_file=${input_location}/celluloidXY/v0.11.7/solution1/segments_${sample}.txt.sorted.bed.gz
sbs_location=${input_location}/cosmicSigNNLS/${sample}_SBS_signatures.txt
SV_vcf_location=${work_dir}/sv.strict/${sample}.sv.cnv.txt

outpath=${work_dir}/hrdetect/

if [ -e "${work_dir}/hrdetect/${sample}.hrdetect.done" ]; then
      echo "file exists, skipping"
else
    
    mkdir ${outpath}
    touch ${work_dir}/hrdetect/${sample}.hrdetect.chk

    Rscript ${script_dir}/sigTools_runthrough.R -t ${sample} --snvFile ${sbs_location} --indelFile ${indel_vcf_location} --SVFile ${SV_vcf_location} --LOHFile ${LOH_seg_file} --outpath ${outpath}/${sample}.hrd.txt --basedir ${script_dir}/

    if [ -e "${outpath}/${sample}.hrd.txt" ]; then
            rm   ${work_dir}/hrdetect/${sample}.hrdetect.chk
            touch ${work_dir}/hrdetect/${sample}.hrdetect.done
    fi
fi


#### pyclone
celluloid_location=/.mounts/labs/PCSI/users/fbeaudry/more.data/${sample}/celluloid/solution/

input_location=/.mounts/labs/PCSI/pipeline/hg38_production/${donor}/${sample}/wgs/bwa/0.7.17

if [ -e "${work_dir}/pyclone/${sample}.pyclone.done" ]; then
      echo "file exists, skipping"
else
    
    outpath=${work_dir}/pyclone/
    mkdir ${outpath}/
    touch ${work_dir}/pyclone/${sample}.pyclone.chk

    if [ ! -e "${celluloid_location}/segments_${sample}.txt.gz" ]; then
        gzip ${celluloid_location}/segments_${sample}.txt
    fi

    python3 ${script_dir}/parse_SNV_freqs.py ${input_location}/final_strelka2_mutect2/${sample}_merged_final.vcf.gz ${outpath}/${sample} ${celluloid_location}/parameters_${sample}.txt ${celluloid_location}/segments_${sample}.txt.gz
    
    if test -f ${input_location}/final_indel/${sample}_merged_final_indels.vcf.gz
    then
    python3  ${script_dir}/parse_indel_freqs.py ${input_location}/final_indel/${sample}_merged_final_indels.vcf.gz ${outpath}/${sample} ${celluloid_location}/parameters_${sample}.txt ${celluloid_location}/segments_${sample}.txt.gz
    else
    python3  ${script_dir}/parse_indel_freqs.py ${input_location}/final_indel.old/${sample}_merged_final_indels.vcf.gz ${outpath}/${sample} ${celluloid_location}/parameters_${sample}.txt ${celluloid_location}/segments_${sample}.txt.gz

    fi

    cp ${outpath}/${sample}.snv.depths.txt ${outpath}/${sample}.all.depths.txt 
    awk '$1 !~ "mutation_id" {print}' ${outpath}/${sample}.dbs.depths.txt >>${outpath}/${sample}.all.depths.txt
    awk '$1 !~ "mutation_id" {print}' ${outpath}/${sample}.indel.depths.txt >>${outpath}/${sample}.all.depths.txt
    awk '$1 !~ "mutation_id" {print}' ${work_dir}/sv.strict/${sample}.sv.depths.txt >>${outpath}/${sample}.all.depths.txt

    pyclone-vi fit -i ${outpath}/${sample}.all.depths.txt -o ${outpath}/${sample}.all.depths.txt.h5 -c 40 -d beta-binomial -r 10

    pyclone-vi write-results-file -i ${outpath}/${sample}.all.depths.txt.h5 -o  ${outpath}/${sample}.all.pyclone.tsv

    Rscript ${script_dir}/plot.pyclone.R --pyclone.path  ${outpath}/${sample}.all.pyclone.tsv --depths.path ${outpath}/${sample}.all.depths.txt --out.path ${outpath}/${sample}


    if [ -e "${outpath}/${sample}.all.pyclone.tsv" ]; then
            rm   ${work_dir}/pyclone/${sample}.pyclone.chk
            touch ${work_dir}/pyclone/${sample}.pyclone.done
    fi
fi


