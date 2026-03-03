#!/bin/bash
#SBATCH -J minimap
#SBATCH -p all
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --array=0
##SBATCH --mem=100G #[M|G|T]
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
REF="/syn2/zhaolian/3.JiLab/reference/PacBio_Targetonly_476bp.mmi"
WORKDIR="/syn2/zhaolian/3.JiLab/results/3.PacBio/1.SMALTout/3.C4007-C4011/01.bam/"
INDIR="/syn2/zhaolian/3.JiLab/results/3.PacBio/1.SMALTout/3.C4007-C4011/02.fastq/"
OUTDIR="/syn2/zhaolian/3.JiLab/results/3.PacBio/1.SMALTout/3.C4007-C4011/03.bam_split/"
mkdir -p $OUTDIR
cd $INDIR

declare -a outname=('C4007-C4011')
describer=${outname[$SLURM_ARRAY_TASK_ID]}
echo ${describer}

#minimap2 -d PacBio_Targetonly_476bp.mmi PacBio_Targetonly_476bp.fasta

minimap2 -t 28 -A 4 -B 12 -O 10,15 -E 2,1 --score-N 0 --end-bonus 10 -a --MD -x map-pb $REF ${describer}.fastq > $OUTDIR/${describer}.sam
samtools view -@ 10 -bS $OUTDIR/${describer}.sam > $OUTDIR/${describer}_unsorted.bam
samtools sort -@ 10 -o $OUTDIR/${describer}.bam -O bam -T ${describer} $OUTDIR/${describer}_unsorted.bam
samtools index $OUTDIR/${describer}.bam
rm $OUTDIR/${describer}.sam
rm $OUTDIR/${describer}_unsorted.bam

#################################
endTime=`date +%Y%m%d-%H:%M:%S`
endTime_s=`date +%s`

sumTime=$[ $endTime_s - $startTime_s ]

echo "$startTime ---> $endTime" "Total:$sumTime seconds"
