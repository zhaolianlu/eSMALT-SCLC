import pandas as pd
import cassiopeia as cas
import os
import cassiopeia
import gzip
import sys
import re
import numpy as np
import sys, threading
sys.setrecursionlimit(10**7) # max depth of recursion
threading.stack_size(2**27)
#########################################################################################
## run main scripts
#########################################################################################



name = sys.argv[1]  #'test_sample'
fastq_path = sys.argv[2]  #'/data/liangzhen/jinhua_jilab_project/data/DNA_Amplicon/T1/staticBC_{}.fastq.gz'
outdir = sys.argv[3]  #'/data2/kantian/LineageTracing/JiLab/01.Target_Q25_Q186_Q192_C4007/01.cassiopeia_out'
print(name)
print(fastq_path)
print(outdir)

input_files = [fastq_path.format('R1'), fastq_path.format('R2')]

print(input_files)

output_directory = outdir# +'/' + name
n_threads = 20
verbose = True
cassiopeia.pp.setup(output_directory, verbose=verbose)

# convert
bam_fp = cas.pp.convert_fastqs_to_unmapped_bam(
	input_files,
	chemistry='10xv3',
	output_directory=output_directory,
	name=name,
	n_threads=n_threads
)
#bam_fp: '/syn2/zhaolian/3.JiLab/Results/1.BarcodeSeq/02.cassiopeia_out/Q192/LvT_unmapped.bam'

# filter_bam
bam_fp = cas.pp.filter_bam(
	bam_fp,
	output_directory=output_directory,
	quality_threshold=10,
	n_threads=n_threads,
)
#bam_fp: '/syn2/zhaolian/3.JiLab/Results/1.BarcodeSeq/02.cassiopeia_out/Q192/LLT_unmapped_filtered.bam'

# error_correct_cellbcs_to_whitelist
bam_fp = cas.pp.error_correct_cellbcs_to_whitelist(
  	bam_fp,
	whitelist='/syn2/zhaolian/3.JiLab/results/1.BarcodeSeq/3M-february-2018.txt',
	output_directory=output_directory,
	n_threads=n_threads,
)
#bam_fp: '/syn2/zhaolian/3.JiLab/Results/1.BarcodeSeq/02.cassiopeia_out/Q192/LLT_unmapped_filtered_corrected.bam'

# collapse
umi_table = cas.pp.collapse_umis(
	bam_fp,
	output_directory=output_directory,
	max_hq_mismatches=3,
	max_indels=2,
	method='likelihood',
	n_threads=n_threads,
)
#bam_fp: '/syn2/zhaolian/3.JiLab/Results/1.BarcodeSeq/02.cassiopeia_out/Q192/LLT_unmapped_filtered_corrected_sorted.bam'

#umi_table = umi_table = pd.read_csv(output_directory+'/'+ name + '_unmapped_filtered_corrected_sorted.collapsed.txt',header=0,sep='\t')

#resolve
umi_table = cas.pp.resolve_umi_sequence(
    umi_table,
    output_directory=output_directory,
    min_umi_per_cell=2,
    min_avg_reads_per_umi=2.0,
    plot=True,
)


umi_table['intBC'] = name 

umi_table = cas.pp.error_correct_umis(
    umi_table,
    max_umi_distance=1,
    allow_allele_conflicts=False,
    n_threads=n_threads,
)

umi_table.to_csv(output_directory + '/'+name+'_umi_table_filtered.csv')

