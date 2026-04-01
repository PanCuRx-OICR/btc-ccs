#!/usr/bin/env bash

# tumour id
#put in argeuments instead
sample="$1"
normal="$2"
donor="$3"

# Check if is provided
if [ -z "$sample" ]; then
    echo "Usage: $0 <sample> <normal> <donor>"
    exit 1
fi
#usage: qsub -P pcsi -l h_vmem=4G,h_rt=3:0:0:0 -cwd -V  -N ${sample}.msi -o ${work_dir}/msi/${sample}/${sample}.msi.out.log -e ${work_dir}/msi/${sample}/${sample}.msi.err.log msisensor.sh ${sample}

##data file setup
#msisensor-pro scan -d /.mounts/labs/PCSI/references/GRCh38/dna/fasta/GRCh38.fa -o /.mounts/labs/PCSI/users/fbeaudry/GRCh38.msi.list

module load msisensorpro/1.2.0

core_dir=/.mounts/labs/PCSI/pipeline/hg38_production
script_dir=/.mounts/labs/PCSI/users/fbeaudry/
work_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/${sample}/msi/
#does it exist, if not create
cd $work_dir

tumor_cram=${core_dir}/${donor}/${sample}/wgs/bwa/0.7.17/collapsed/${sample}.cram
normal_cram=${core_dir}/${donor}/${normal}/wgs/bwa/0.7.17/collapsed/${normal}.cram

tumor_bam=${core_dir}/${donor}/${sample}/wgs/bwa/0.7.17/collapsed/${sample}.bam
normal_bam=${core_dir}/${donor}/${normal}/wgs/bwa/0.7.17/collapsed/${normal}.bam

if [ -e "${work_dir}/${sample}.msi.done" ]; then
      echo "file exists, skipping"
else
      touch ${work_dir}/${sample}.msi.chk

      for boot in {1..20}
      do

      shuf -n 500 ${script_dir}/GRCh38.msi.list >${work_dir}/rep.list

      awk '$2 !~ "location" {print}' ${work_dir}/rep.list | \
            sort -k1,1 -k2,2n >${work_dir}/rep.list.sorted

      if test -f $tumor_cram
      then

      msisensor-pro msi \
            -d ${work_dir}/rep.list.sorted \
            -n ${normal_cram} -t ${tumor_cram} \
            -g /.mounts/labs/PCSI/references/GRCh38/dna/fasta/GRCh38.fa \
            -o ${work_dir}/${sample}.msi | tee -a ${work_dir}/${sample}.msi.log

      

      awk -v boot="${boot}" '$1 !~ "Total_Number_of_Sites" {print boot"\t"$1"\t"$2"\t"$3}' ${work_dir}/${sample}.msi >>${work_dir}/${sample}.msi.booted

      #check if error code is not zero

      else

      #add more logging statements

      msisensor-pro msi \
            -d ${work_dir}/rep.list.sorted \
            -n ${normal_bam} -t ${tumor_bam} \
            -g /.mounts/labs/PCSI/references/GRCh38/dna/fasta/GRCh38.fa \
            -o ${work_dir}/${sample}.msi

            #add tee thing

      awk -v boot="${boot}" '$1 !~ "Total_Number_of_Sites" {print boot"\t"$1"\t"$2"\t"$3}' ${work_dir}/${sample}.msi >>${work_dir}/${sample}.msi.booted

      fi

      done

      if [ -e "${work_dir}/${sample}.msi.booted" ]; then
            rm   ${work_dir}/${sample}.msi.chk
            touch ${work_dir}/${sample}.msi.done
      fi
fi