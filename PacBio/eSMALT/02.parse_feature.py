import multiprocessing as mp
import os
import sys
from dataclasses import dataclass
from typing import Tuple

import edlib
import pysam
import utils
 
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
@dataclass
class Seq:
    seq: str = ""
    qual: str = ""

@dataclass
class MatchedSeq:
    #P5: Seq
    cellBC: Seq
    umi: Seq
    bio: Seq
    intBC: Seq
    i7: Seq
    i5: Seq

def lookup_pos(pos, span_dict=None):
    """
    3KYC primers
    [10:19] : P5
    [34:49] : cell barcode =>
    [50:61] : UMI  #UMI [X, X+3]
    [355:830] : bio =>
    [848:874] : static barcode <=
    """

    span_dict = {
        # range(9, 19): "P5",
        # range(33, 49): "cellBC",
        # range(49, 61+3): "umi",
        # range(354, 830): "bio",
        # range(847, 874): "intBC",
        range(914, 930): "cellBC",
        range(902, 914): "umi",
        range(133, 609): "bio",
        range(89, 116): "intBC",
        range(24, 34): "i7",
        range(963, 973): "i5",
    }
    for span in span_dict:
        if pos in span:
            return span_dict[span]

def parse_read(read):
    """parse bioseq to fastq."""
    #  dele = "-"
    dele = ""

    # check if there is MD tag
    # check whether query sequence is exist
    if not read.has_tag("MD") or (query_seq := read.query_sequence) is None:
        return None

    query_qual = read.query_qualities
    matched = MatchedSeq(
        #P5=Seq(seq="", qual=""),
        cellBC=Seq(seq="", qual=""),
        umi=Seq(seq="", qual=""),
        bio=Seq(seq="", qual=""),
        intBC=Seq(seq="", qual=""),
        i7=Seq(seq="", qual=""),
        i5=Seq(seq="", qual=""),
    )
    ref_pos_pointer = -1
    for query_pos, ref_pos, ref_base in read.get_aligned_pairs(with_seq=True):

        if ref_pos is not None:
            ref_pos_pointer = ref_pos
        feature_type = lookup_pos(ref_pos_pointer)
        if feature_type in ["cellBC"] and ref_base is not None:
            getattr(matched, feature_type).seq += (query_seq[query_pos] if query_pos is not None else dele)
        elif feature_type in ["umi","bio","intBC","i7","i5"]:
            getattr(matched, feature_type).seq += (query_seq[query_pos] if query_pos is not None else dele)
            getattr(matched, feature_type).qual += (
                chr(query_qual[query_pos] + 33) ## 33??
                if query_pos is not None
                else dele
            )

    if (
        (len(matched.bio.seq) > 300)
#         and check_library_barcode(matched.lib5.seq, matched.lib3.seq)
        and (matched.cellBC.seq is not None) and (matched.intBC.seq is not None)
    ):
        
        return (
            matched.cellBC.seq,
            matched.umi.seq,
            matched.umi.qual,
            matched.intBC.seq,
            matched.intBC.qual,
            matched.bio.seq,
            matched.bio.qual,
            matched.i7.seq,
            matched.i7.qual,
            matched.i5.seq,
            matched.i5.qual,
        )
    return None




#########################################################################################
## 1. run main scripts
#########################################################################################
if __name__ == "__main__":
    """
    2nd speed-limit step(2/2): 150 G data costs about 8 hours
    """
    SAMPLE_NAME = sys.argv[1]      # SAMPLE_NAME = "LLT-1"
    INDIR=sys.argv[2]
    OUTDIR=sys.argv[3]

    
    INPUT_BAM_FILE = INDIR+SAMPLE_NAME+".bam"
    OUT_FASTQ=OUTDIR+SAMPLE_NAME+".fastq"
    OUT_TEXT=OUTDIR+SAMPLE_NAME+".txt"
    if not os.path.exists(OUTDIR):
        os.mkdir(OUTDIR)

    with pysam.AlignmentFile(
        INPUT_BAM_FILE, "rb", check_sq=False
    ) as infile, open(OUT_FASTQ,"w") as outq, open(OUT_TEXT,"w") as outt :
#         f={}
        for read in infile.fetch(until_eof=True):
            result = parse_read(read)
            if result:
                outq.write(
                    f"@{read.qname}\n{result[5]}\n+\n{result[6]}\n"
                )
                outt.write(
                    f"@{read.qname}\t{result[0]}\t{result[1]}\t{result[2]}\t{result[3]}\t{result[4]}\t{result[7]}\t{result[8]}\t{result[9]}\t{result[10]}\n"
                )
