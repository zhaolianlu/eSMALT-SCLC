#!/bin/bash
#SBATCH -J bwaC4007
#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 5
#SBATCH --array=0-7
#SBATCH --mem=10G #[M|G|T]
#SBATCH -t 1000:00:00 #[days-hh:mm:ss]
#SBATCH -o ./temp/%x-%j.log #%x: job_name, %j: job_id
#SBATCH -e ./temp/%x-%j.err #leave out -e flag to join error/log files 
##SBATCH --mail-user=zl.lu@siat.ac.cn
##SBATCH --mail-type=BEGIN,END #common valid events: NONE, BEGIN, END, FAIL, REQUEUE, ALL

startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`
################################

FASTQDIR='/syn2/zhaolian/3.JiLab/results/1.BarcodeSeq/01.fastqMerged/'
OUTDIR='/syn2/zhaolian/3.JiLab/results/1.BarcodeSeq/02.bwa_out/C4007/'
REF='/syn2/zhaolian/3.JiLab/reference/PacBio_Targetonly_476bp.fasta'
mkdir -p ${OUTDIR}

declare -a array=('C4007-LLT-1_merged' 'C4007-LLT-2_merged' 'C4007-Abdo_merged' 'C4007-Chest_merged' 'C4007-LLT-1_mergednew' 'C4007-LLT-2_mergednew' 'C4007-Abdo_mergednew' 'C4007-Chest_mergednew')
declare -a labels=('LLT-1-old' 'LLT-2-old' 'Abdo-old' 'Chest-old' 'LLT-1' 'LLT-2' 'Abdo' 'Chest')
inname=${array[$SLURM_ARRAY_TASK_ID]}
outname=${labels[$SLURM_ARRAY_TASK_ID]}
echo $inname


cd $FASTQDIR

## 0.make index, run once
#bwa index -a bwtsw CB_UMI_polyT_target476_SB_P862bp.fasta

## 1.bwa
bwa mem -M -Y -t 5 ${REF} ${inname}_R1.fastq.gz ${inname}_R2.fastq.gz > ${OUTDIR}/${outname}.sam

## 2.samtools
cd $OUTDIR
samtools view -bS ${outname}.sam > ${outname}.uns.bam
samtools sort -@ 5 -n ${outname}.uns.bam > ${outname}.bam
#samtools index ${outname}.bam
rm ${outname}.uns.bam
samtools flagstat ${outname}.bam > ${outname}_stats.txt

## 3. extract mapped reads

sambamba view -f bam -F "not (unmapped and mate_is_unmapped)" ${outname}.bam > ${outname}_sambamba.bam
sambamba view -f bam -F "not secondary_alignment" ${outname}_sambamba.bam > ${outname}_sambamba_mapped.bam
rm ${outname}_sambamba.bam

## 4. convert to fastq
samtools bam2fq -n ${outname}_sambamba_mapped.bam -1 ${outname}_mapped_R1.fastq -2 ${outname}_mapped_R2.fastq
gzip ${outname}_mapped_R1.fastq
gzip ${outname}_mapped_R2.fastq


#################################
endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`

sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"




