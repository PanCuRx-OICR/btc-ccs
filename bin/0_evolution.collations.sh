

conda activate pyclone-vi
data_dir=/.mounts/labs/PCSI/users/fbeaudry/
work_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/
sample_file=more.sample.list.txt

mysamples=$( awk '{ print $2 }' $data_dir/${sample_file} | tr '\n' ' ' )

for tumor in ${mysamples[@]}
do

echo ${tumor}

mkdir -p ${work_dir}/${tumor}/sex/

donor=$(awk -v tumor=${tumor} '$2==tumor {print $1}' $data_dir/${sample_file} )
normal_sample=$(awk -v tumor=${tumor} '$2==tumor {print $3}' $data_dir/${sample_file} )

rm ${work_dir}/${tumor}/sex/${normal_sample}.out.log ${work_dir}/${tumor}/sex/${normal_sample}.err.log 
qsub -P pcsi -l h_vmem=6G,h_rt=0:12:0:0 -cwd -V  -N ${tumor}.sex -o ${work_dir}/${tumor}/sex/${normal_sample}.out.log -e ${work_dir}/${tumor}/sex/${normal_sample}.err.log /.mounts/labs/PCSI/users/fbeaudry/more.data/btc.scripts/xy_coverage.sh ${normal_sample} ${donor} ${tumor} ${work_dir}/

done 


#### celluloid

conda activate pyclone-vi

script_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/btc.scripts/
data_dir=/.mounts/labs/PCSI/users/fbeaudry/
work_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/
sample_file=more.sample.list.txt

mysamples=$( awk '{ print $2 }' $data_dir/${sample_file} | tr '\n' ' ' )

for tumor in ${mysamples[@]}
do

donor=$(awk -v tumor=${tumor} '$2==tumor {print $1}' $data_dir/${sample_file} )
normal_sample=$(awk -v tumor=${tumor} '$2==tumor {print $3}' $data_dir/${sample_file} )

work_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/
cell_dir=${work_dir}/${tumor}/celluloid/
mkdir -p ${cell_dir}
cd  ${cell_dir}

rm ${cell_dir}/${tumor}.celluloid.out.log ${cell_dir}/${tumor}.celluloid.err.log
qsub -P pcsi -l h_vmem=6G,h_rt=2:0:0:0 -cwd -V  -N ${tumor}.cell -o ${cell_dir}/${tumor}.celluloid.out.log -e ${cell_dir}/${tumor}.celluloid.err.log ${script_dir}/launch_celluloid.sh ${tumor} ${normal_sample} ${donor} ${work_dir}
done

### add text file with final celluloid solution

for tumor in ${mysamples[@]}
do

#readlink /.mounts/labs/PCSI/users/fbeaudry/more.data/${tumor}/celluloid/solution >/.mounts/labs/PCSI/users/fbeaudry/more.data/${tumor}/celluloid/solution.txt
cat /.mounts/labs/PCSI/users/fbeaudry/more.data/${tumor}/celluloid/solution.txt

done



#### the rest
conda activate pyclone-vi

work_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/
cd ${work_dir}

data_dir=/.mounts/labs/PCSI/users/fbeaudry/
sample_file=more.sample.list.txt

mysamples=$( awk '{ print $2 }' $data_dir/${sample_file} | tr '\n' ' ' )

for tumor in ${mysamples[@]}
do

echo ${tumor}

mkdir -p ${work_dir}/${tumor}/

donor=$(awk -v tumor=${tumor} '$2==tumor {print $1}' $data_dir/${sample_file} )
normal_sample=$(awk -v tumor=${tumor} '$2==tumor {print $3}' $data_dir/${sample_file} )

rm ${work_dir}/${tumor}/${tumor}.out.log ${work_dir}/${tumor}/${tumor}.err.log

#rm ${work_dir}/${tumor}/sv.strict/${tumor}.sv_strict.done
qsub -P pcsi -l h_vmem=6G,h_rt=2:0:0:0 -cwd -V  -N ${tumor}.extra -o ${work_dir}/${tumor}/${tumor}.out.log -e ${work_dir}/${tumor}/${tumor}.err.log /.mounts/labs/PCSI/users/fbeaudry/more.data/btc.scripts/fe.extra.pipe.sh ${tumor} ${normal_sample} ${donor}

#### launch MSI ####
mkdir ${work_dir}/${tumor}/msi/
rm ${work_dir}/${tumor}/msi/${tumor}.msi.out.log ${work_dir}/${tumor}/msi/${tumor}.msi.err.log
qsub -P pcsi -l h_vmem=4G,h_rt=0:20:0:0 -cwd -V  -N ${tumor}.msi -o ${work_dir}/${tumor}/msi/${tumor}.msi.out.log -e ${work_dir}/${tumor}/msi/${tumor}.msi.err.log /.mounts/labs/PCSI/users/fbeaudry/more.data/btc.scripts/msisensor.sh ${tumor} ${normal_sample} ${donor}

done 







#### sigprofiler ####

conda activate /.mounts/labs/PCSI/users/fbeaudry/py_envs/sigprofiler

out_root=/.mounts/labs/PCSI/users/fbeaudry/more.data/
script_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/btc.scripts/

data_dir=/.mounts/labs/PCSI/users/fbeaudry/
sample_file=sample.list.txt

mysamples=$( awk '{ print $2 }' $data_dir/${sample_file} | tr '\n' ' ' )

for tumour_id in ${mysamples[@]}
do

echo $tumour_id
donor=$(awk -v tumor=${tumour_id} '$2==tumor {print $1}' $data_dir/${sample_file} )

work_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/${tumour_id}/

mkdir -p ${work_dir}/sigprofiler/

rm ${work_dir}/sigprofiler/${tumour_id}.sig.out.log ${work_dir}/sigprofiler/${tumour_id}.sig.err.log 

qsub -P pcsi -l h_vmem=6G,h_rt=2:0:0:0 -cwd -V  -N ${tumour_id}.sig -o ${work_dir}/sigprofiler/${tumour_id}.sig.out.log -e ${work_dir}/sigprofiler/${tumour_id}.sig.err.log /.mounts/labs/PCSI/users/fbeaudry/more.data/btc.scripts/sigprofiler.sh ${tumour_id} ${donor}

done


#### collate signatures

data_dir=/.mounts/labs/PCSI/users/fbeaudry/

mysamples=$( awk '{ print $2 }' $data_dir/sample.list.txt | tr '\n' ' ' )

sample=BTC_0021_Gb_P_526
donor=$(awk -v tumor=${sample} '$2==tumor {print $1}' $data_dir/sample.list.txt )

awk '$1 ~ "Samples" {print}' /.mounts/labs/PCSI/users/fbeaudry/more.data/${sample}/sigprofiler/96/Assignment_Solution/Activities/Assignment_Solution_Activities.txt >/.mounts/labs/PCSI/users/fbeaudry/more.data/SBS.sig.txt
awk '$1 ~ "Samples" {print}' /.mounts/labs/PCSI/users/fbeaudry/more.data/${sample}/sigprofiler/DINUC/Assignment_Solution/Activities/Assignment_Solution_Activities.txt >/.mounts/labs/PCSI/users/fbeaudry/more.data/DBS.sig.txt
awk '$1 ~ "Samples" {print}' /.mounts/labs/PCSI/users/fbeaudry/more.data/${sample}/sigprofiler/ID/Assignment_Solution/Activities/Assignment_Solution_Activities.txt >/.mounts/labs/PCSI/users/fbeaudry/more.data/indel.sig.txt

for sample in ${mysamples[@]}
do

echo ${sample}

donor=$(awk -v tumor=${sample} '$2==tumor {print $1}' $data_dir/sample.list.txt )

awk '$1 !~ "Samples" {print}' /.mounts/labs/PCSI/users/fbeaudry/more.data/${sample}/sigprofiler/96/Assignment_Solution/Activities/Assignment_Solution_Activities.txt >>/.mounts/labs/PCSI/users/fbeaudry/more.data/SBS.sig.txt

awk '$1 !~ "Samples" {print}' /.mounts/labs/PCSI/users/fbeaudry/more.data/${sample}/sigprofiler/DINUC/Assignment_Solution/Activities/Assignment_Solution_Activities.txt >>/.mounts/labs/PCSI/users/fbeaudry/more.data/DBS.sig.txt

awk '$1 !~ "Samples" {print}' /.mounts/labs/PCSI/users/fbeaudry/more.data/${sample}/sigprofiler/ID/Assignment_Solution/Activities/Assignment_Solution_Activities.txt >>/.mounts/labs/PCSI/users/fbeaudry/more.data/indel.sig.txt


done

#### mutationtimer

data_dir=/.mounts/labs/PCSI/users/fbeaudry/
sample_file=sample.list.txt

mysamples=$( awk '{ print $2 }' $data_dir/${sample_file} | tr '\n' ' ' )

for sample in ${mysamples[@]}

do
work_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/${sample}/


mkdir ${work_dir}/mutationtimer/
cat ${work_dir}/pyclone/${sample}.indel.mutationTimeR.vcf ${work_dir}/sv.strict/${sample}.sv.mutationTimeR.headless.vcf >${work_dir}/mutationtimer/${sample}.sv.mutationTimeR.vcf
awk '$7 ~ "PASS" {print}' ${work_dir}/pyclone/${sample}.dbs.mutationTimeR.vcf >>${work_dir}/mutationtimer/${sample}.sv.mutationTimeR.vcf

cp /.mounts/labs/PCSI/users/fbeaudry/more.data/${sample}/celluloid/solution/*.txt ${work_dir}/mutationtimer/
cp /.mounts/labs/PCSI/users/fbeaudry/more.data/${sample}/celluloid/solution/*.txt.gz ${work_dir}/mutationtimer/
## rest runs locally in 4.2_make_evolution_file.R

done


#### ampliconArchitect ####

conda activate ampsuite

script_dir=/.mounts/labs/PCSI/users/fbeaudry/
data_dir=/.mounts/labs/PCSI/users/fbeaudry/

sample_list=$data_dir/sample.list.txt

mysamples=$( awk '{ print $2 }' $sample_list | tr '\n' ' ' )

for tumor in ${mysamples[@]}
do

echo ${tumor}

work_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/${tumor}
outpath=${work_dir}/amplicon/
mkdir ${outpath}

donor=$(awk -v tumor=${tumor} '$2==tumor {print $1}' $sample_list )

rm ${outpath}/${tumor}.amplicon.out.log ${outpath}/${tumor}.amplicon.err.log 

qsub -P pcsi -l h_vmem=8G,h_rt=18:0:0:0  -pe smp 12 -cwd -V  -N ${tumor}.amplicon -o ${outpath}/${tumor}.amplicon.out.log -e ${outpath}/${tumor}.amplicon.err.log ${script_dir}/more.data/btc.scripts/call_aa.sh ${tumor} ${donor}

done

#### run gistic ####

#cd /.mounts/labs/PCSI/users/fbeaudry/gistic/MCR_Installer
# ./install -mode silent -agreeToLicense yes -destinationFolder /.mounts/labs/PCSI/users/fbeaudry/gistic/MATLAB_Compiler_Runtime



script_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/celluloid/
data_dir=/.mounts/labs/PCSI/users/fbeaudry/
refgenefile=/.mounts/labs/PCSI/users/fbeaudry/gistic/refgenefiles/hg38.UCSC.add_miR.160920.refgene.mat

awk '{print $2}' $data_dir/sample.list.txt >$data_dir/all.list.txt
 
for set in  CMSA CMSB all
do

cd $data_dir

data_dir=/.mounts/labs/PCSI/users/fbeaudry/
sample_file=$data_dir/${set}.list.txt

segfile=/.mounts/labs/PCSI/users/fbeaudry/more.data/${set}.seg.txt
out_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/gistic/${set}/
mkdir $out_dir

mysamples=$( awk '{ print $1 }' $sample_file | tr '\n' ' ' )

tumor=BTC_9007_Lv_P_526
work_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/${tumor}/
zcat ${work_dir}/celluloid/solution/segments_${tumor}.txt.gz | awk '$1 ~ "chrom" {print}' >/.mounts/labs/PCSI/users/fbeaudry/more.data/${set}.cnv.txt

rm ${segfile} 
for tumor in ${mysamples[@]}
do
echo ${tumor}

work_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/${tumor}/

awk '$7 ~ "FALSE" {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' ${work_dir}/celluloid/solution/${tumor}.classic.seg  >>${segfile}
zcat ${work_dir}/celluloid/solution/segments_${tumor}.txt.gz | awk '$1 !~ "chrom" {print}' >>/.mounts/labs/PCSI/users/fbeaudry/more.data/${set}.cnv.txt

done 

#from head node (trouble on compute)
cd /.mounts/labs/PCSI/users/fbeaudry/gistic/
/.mounts/labs/PCSI/users/fbeaudry/gistic/gistic2 -b $out_dir -seg $segfile -refgene $refgenefile -genegistic 1 -broad 0 -armpeel 0 -savegene 1 -conf 0.9
done

#### celullarity stats ####

data_dir=/.mounts/labs/PCSI/users/fbeaudry/

mysamples=$( awk -F'\t' '{ print $3 }' $data_dir/every.sample.list.txt | tr '\n' ' ' )

rm ${data_dir}/btc.coverage.txt ${data_dir}/btc.cellularity.curated.txt ${data_dir}/btc.cellularity.pipeline.txt ${data_dir}/btc.rna.coverage.txt
for sample in ${mysamples[@]}
do

echo ${sample}

do_nor=$(awk -F'\t' -v tumor=${sample} '$3==tumor {print $1}' $data_dir/every.sample.list.txt )
donor=$(echo ${do_nor} | tr -d _)

sed -n '8p' /.mounts/labs/PCSI/pipeline/hg38_production/${donor}/${sample}/wgs/bwa/0.7.17/coverage/${sample}_coverage_collapsed.metrics | awk -v sample=${sample} '{print sample, "T", $2, $4}' >>${data_dir}/btc.coverage.txt

cat /.mounts/labs/PCSI/users/fbeaudry/more.data/${sample}/celluloid/solution/parameters_${sample}.txt | awk -v sample=${sample} '$1 !~ "value" {print sample, $4, $5}' >>${data_dir}/btc.cellularity.curated.txt
cat /.mounts/labs/PCSI/pipeline/hg38_production/${donor}/${sample}/wgs/bwa/0.7.17/celluloidXY/v0.11.7/solution/parameters_${sample}.txt | awk -v sample=${sample} '$1 !~ "value" {print sample, $4, $5}' >>${data_dir}/btc.cellularity.pipeline.txt

cat /.mounts/labs/PCSI/pipeline/hg38_production/${donor}/${sample}/rna/star/2.7.4a/collapsed/Log.final.out | awk -F'|' -v tumor=${sample} '$1 ~ "Uniquely mapped reads number" {print tumor"\t"$2}' >>${data_dir}/btc.rna.coverage.txt

done

## normal
data_dir=/.mounts/labs/PCSI/users/fbeaudry/
sample_file=sample.list.txt

mysamples=$( awk '{ print $2 }' $data_dir/${sample_file} | tr '\n' ' ' )

for tumor in ${mysamples[@]}
do

donor=$(awk -v tumor=${tumor} '$2==tumor {print $1}' $data_dir/${sample_file} )
normal_sample=$(awk -v tumor=${tumor} '$2==tumor {print $3}' $data_dir/${sample_file} )

echo ${normal_sample}
sed -n '8p' /.mounts/labs/PCSI/pipeline/hg38_production/${donor}/${normal_sample}/wgs/bwa/0.7.17/coverage/${normal_sample}_coverage_collapsed.metrics | awk -v sample=${normal_sample} '{print sample, "N", $2, $4}' >>${data_dir}/btc.coverage.txt

done




#### DDR genes ####

data_dir=/.mounts/labs/PCSI/users/fbeaudry/
work_dir=/.mounts/labs/PCSI/users/fbeaudry
core_dir=/.mounts/labs/PCSI/pipeline/hg38_production

mysamples=$( awk '{ print $2 }' $data_dir/sample.list.txt | tr '\n' ' ' )

donor=BTC0027
sample=BTC_0027_Lv_M_526

awk '$1 ~ "donor" {print}' ${core_dir}/${donor}/${sample}/wgs/bwa/0.7.17/results/${sample}.variants.csv |   awk -F','  '{print}'  > $work_dir/variants.ddr.csv

for sample in ${mysamples[@]}
do
echo $sample
donor=$(awk -v sample=${sample} '$2==sample {print $1}' $data_dir/sample.list.txt )

while read this_gene
do

echo $this_gene


if test -f ${core_dir}/${donor}/${sample}/wgs/bwa/0.7.17/results/${sample}.variants.csv
then
awk '$1 !~ "donor" {print}' ${core_dir}/${donor}/${sample}/wgs/bwa/0.7.17/results/${sample}.variants.csv | grep $this_gene |  awk -F',' -v this_gene=$this_gene '$10 == this_gene {print}'  >>  $work_dir/variants.ddr.csv
else
awk '$1 !~ "donor" {print}' ${core_dir}/${donor}/${sample}/wgs/bwa/0.7.17/results.old/${sample}.variants.csv | grep $this_gene |  awk -F',' -v this_gene=$this_gene '$10 == this_gene {print}'  >>  $work_dir/variants.ddr.csv
fi

done < ${data_dir}/ddr.genes
done



#### fusions ####

core_dir=/.mounts/labs/PCSI/pipeline/hg38_production
work_dir=/.mounts/labs/PCSI/users/fbeaudry
data_dir=/.mounts/labs/PCSI/users/fbeaudry/

cd $work_dir

#mysamples=$( awk '{ print $2 }' $work_dir/sample.list.txt | tr '\n' ' ' )
mysamples=$( awk -F'\t' '{ print $3 }' $data_dir/every.sample.list.txt | tr '\n' ' ' )

rm $work_dir/BTC.star.fusions.txt $work_dir/BTC.mavis.fusions.txt

for sample in ${mysamples[@]}
do
do_nor=$(awk -F'\t' -v tumor=${sample} '$3==tumor {print $1}' $data_dir/every.sample.list.txt )
donor=$(echo ${do_nor} | tr -d _)

echo $donor
awk -v donor=$donor -v sample=$sample '$1 !~ "#FusionName" {print donor"\t"sample"\t"$1}' ${core_dir}/${donor}/${sample}/rna/star/2.7.4a/starfusion/1.9.0/star-fusion.fusion_predictions.tsv >>$work_dir/BTC.star.fusions.txt
awk -F '\t' -v donor=$donor -v sample=$sample '$1 !~ "#tracking_id" {print donor"\t"sample"\t"$5"\t"$10"\t"$11"\t"$12}' ${core_dir}/${donor}/${sample}/rna/star/2.7.4a/starfusion/1.9.0/mavis/summary/mavis_summary_all_${sample}.tab >>$work_dir/BTC.mavis.fusions.txt

done 


#### LOY #####

data_dir=/.mounts/labs/PCSI/users/fbeaudry/
work_dir=/.mounts/labs/PCSI/users/fbeaudry/more.data/

mysamples=$( awk '{ print $2 }' $data_dir/sample.list.txt | tr '\n' ' ' )

rm ${work_dir}/Y_segments.txt ${work_dir}/Y_rna.txt rm ${work_dir}/sex.txt
for tumor in ${mysamples[@]}
do
echo ${tumor}

donor=$(awk -v tumor=${tumor} '$2==tumor {print $1}' $data_dir/sample.list.txt  )
normal=$(awk -v tumor=${tumor} '$2==tumor {print $3}' $data_dir/sample.list.txt  )

echo ${donor}

#celluloid

awk -v tumour_id=$tumor '{print tumour_id"\t"$1}' /.mounts/labs/PCSI/users/fbeaudry/more.data/${tumor}/sex/${normal}.sex.txt >>${work_dir}/sex.txt

#segments
zcat /.mounts/labs/PCSI/users/fbeaudry/more.data/${tumor}/celluloid/solution/segments_${tumor}.txt.gz | awk '$1 ~ "chr23" {print }'  >>${work_dir}/Y_segments.txt
zcat /.mounts/labs/PCSI/users/fbeaudry/more.data/${tumor}/celluloid/solution/segments_${tumor}.txt.gz | awk '$1 ~ "chr24" {print }' >>${work_dir}/Y_segments.txt

#rna
awk -v tumour_id=$tumor '$3 ~ "chrX" {print tumour_id"\t"$2"\t"$7"\t"$8}' /.mounts/labs/PCSI/pipeline/hg38_production/${donor}/${tumor}/rna/star/2.7.4a/stringtie/2.0.6/${tumor}_stringtie_abundance.txt >>${work_dir}/Y_rna.txt
awk -v tumour_id=$tumor '$3 ~ "chrY" {print tumour_id"\t"$2"\t"$7"\t"$8}' /.mounts/labs/PCSI/pipeline/hg38_production/${donor}/${tumor}/rna/star/2.7.4a/stringtie/2.0.6/${tumor}_stringtie_abundance.txt >>${work_dir}/Y_rna.txt


done


##### RNA NMF ####

conda activate /.mounts/labs/PCSI/users/fbeaudry/nmf_env

cd /.mounts/labs/PCSI/users/fbeaudry/nmf_analysis

script_dir=/.mounts/labs/PCSI/users/fbeaudry/nmf_analysis/

rm nmf.out.log.txt nmf.out.err.txt
rm *.rds *.csv
qsub -P pcsi -l h_vmem=60G,h_rt=4:0:0:0 -cwd -V  -N nmf -o nmf.out.log.txt -e nmf.out.err.txt  ${script_dir}/launch_cluster.nmf.sh 

script_dir=/.mounts/labs/PCSI/users/fbeaudry/nmf_analysis/cmsa/
qsub -P pcsi -l h_vmem=60G,h_rt=4:0:0:0 -cwd -V  -N nmf -o ${script_dir}/nmf.out.log.txt -e ${script_dir}/nmf.out.err.txt  ${script_dir}/launch_cluster.nmf.cmsa.sh 

