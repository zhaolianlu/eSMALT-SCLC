#!/bin/bash
#SBATCH -J 4007
#SBATCH -p fat
#SBATCH -N 1
#SBATCH -n 20
##SBATCH --array=0
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
REF="/syn2/zhaolian/2.eSMALT/reference/scAmplicon_16-9-Fi_P5R_1510bp.mmi"
WORKDIR="/syn2/zhaolian/3.JiLab/results/5.Steroseq/01.saw_count"
mkdir -p ${WORKDIR}
cd ${WORKDIR}

#####################################################################################################
# 1. ccs, by pbccs

/syn2/zhaolian/3.JiLab/results/5.Steroseq/SAW/saw-8.1.3/bin/saw count \
    --id=4007 \
    --sn=A05099C413 \
    --omics=transcriptomics \
    --kit-version="Stereo-seq T FF V1.3" \
    --sequencing-type="PE75_50+100" \
    --chip-mask=/syn2/zhaolian/3.JiLab/data/StereoSeq/Rawdata/C4007_PT/A05099C413.barcodeToPos.h5 \
    --organism=mouse \
    --tissue=lung \
    --fastqs=/syn2/zhaolian/3.JiLab/data/StereoSeq/Rawdata/C4007_PT/ \
    --reference=/syn2/zhaolian/3.JiLab/results/5.Steroseq/reference/transcriptome \
    --image-tar=/syn2/zhaolian/3.JiLab/data/StereoSeq/Rawdata/C4007_PT/A05099C413_SC_20250524_111750_4.1.0.tar.gz




#################################
endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`

sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"
