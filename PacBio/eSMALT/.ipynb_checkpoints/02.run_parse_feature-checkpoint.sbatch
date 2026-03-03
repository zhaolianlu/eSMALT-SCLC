#!/bin/bash
#SBATCH -J 02
#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 2
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

declare -a INNAME=('C4007-C4011')
SAMPLE_NAME=${INNAME[$SLURM_ARRAY_TASK_ID]}
INDIR="/syn2/zhaolian/3.JiLab/results/3.PacBio/1.SMALTout/3.C4007-C4011/01.bam/"
OUTDIR="/syn2/zhaolian/3.JiLab/results/3.PacBio/1.SMALTout/3.C4007-C4011/02.fastq/"
mkdir -p ${OUTDIR}
echo ${SAMPLE_NAME}

python 02.parse_feature.py $SAMPLE_NAME $INDIR $OUTDIR

#################################
endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`

sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"
