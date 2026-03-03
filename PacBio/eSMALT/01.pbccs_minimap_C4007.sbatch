#!/bin/bash
#SBATCH -J ccsC4007
#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --array=0
#SBATCH --mem=10G #[M|G|T]
#SBATCH -t 1000:00:00 #[days-hh:mm:ss]
#SBATCH -o ./temp/%x-%j.log #%x: job_name, %j: job_id
#SBATCH -e ./temp/%x-%j.err #leave out -e flag to join error/log files 
##SBATCH --mail-user=zl.lu@siat.ac.cn
##SBATCH --mail-type=BEGIN,END #common valid events: NONE, BEGIN, END, FAIL, REQUEUE, ALL

startTime=`date +%Y%m%d-%H:%M:%S`
startTime_s=`date +%s`
################################

##1st speed-limit step(1/2): 150 G data costs about 24 hours
#####################################################################################################
REF="/syn2/zhaolian/3.JiLab/reference/scAmplicon-SMALT.mmi"
WORKDIR="/syn2/zhaolian/3.JiLab/results/3.PacBio/1.SMALTout/3.C4007-C4011/01.bam"
mkdir -p ${WORKDIR}

declare -a inbams=('/syn2/zhaolian/3.JiLab/data/4.C4007/PacBio/batch3_20240803/m84130_240730_180513_s2.hifi_reads.bc2095.bam')
declare -a outname=('C4007-C4011')
INBAM=${inbams[$SLURM_ARRAY_TASK_ID]}
SAMPLE=${outname[$SLURM_ARRAY_TASK_ID]} 
echo $SAMPLE
echo $INBAM

#####################################################################################################
# 1. ccs, by pbccs
mkdir -p $WORKDIR
cd $WORKDIR
#ccs $INBAM ${SAMPLE}_ccs.bam --min-length=300 --num-threads=20 --by-strand --min-passes=3 --min-rq=0.99 --report-file=${SAMPLE}_ccs.report
samtools fastq -@ 20 $INBAM > ${SAMPLE}_ccs.fastq
rm ${SAMPLE}_ccs.bam
# 2. map, by minimap2
#minimap2 -d scAmplicon-SMALT.mmi scAmplicon-SMALT.fasta
minimap2 -t 20 -A 4 -B 12 -O 10,15 -E 2,1 --score-N 0 --end-bonus 10 -a --MD -x map-pb $REF ${SAMPLE}_ccs.fastq> ${SAMPLE}_aln.sam
samtools view -bS ${SAMPLE}_aln.sam > ${SAMPLE}_unsorted.bam
samtools sort -@ 10 -o ${SAMPLE}.bam -O bam -T ${SAMPLE} ${SAMPLE}_unsorted.bam
samtools index ${SAMPLE}.bam
rm ${SAMPLE}_unsorted.bam
rm ${SAMPLE}_aln.sam

#################################
endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`

sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"
