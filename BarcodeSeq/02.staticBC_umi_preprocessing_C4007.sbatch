#!/bin/bash
#SBATCH -J 02_C4007
#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --array=0-7
#SBATCH --mem=5G #[M|G|T]
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

MOUSE="C4007"
declare -a INNAME=('LLT-1-old' 'LLT-2-old' 'Abdo-old' 'Chest-old' 'LLT-1' 'LLT-2' 'Abdo' 'Chest')

SAMPLE=${INNAME[$SLURM_ARRAY_TASK_ID]}
FQ="/syn2/zhaolian/3.JiLab/results/1.BarcodeSeq/02.bwa_out/${MOUSE}/${SAMPLE}_mapped_{}.fastq.gz"
OUTDIR="/syn2/zhaolian/3.JiLab/results/1.BarcodeSeq/03.cassiopeia_out/${MOUSE}/"
mkdir -p ${OUTDIR}
echo $SAMPLE
echo $FQ
echo $OUTDIR
python 02.staticBC_umi_preprocessing.py $SAMPLE $FQ $OUTDIR

#################################
endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`

sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"
