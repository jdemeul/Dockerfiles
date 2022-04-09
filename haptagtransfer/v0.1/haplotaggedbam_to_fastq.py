#!/usr/bin/env python

import pysam
import argparse
import subprocess

parser = argparse.ArgumentParser(description='make fastq from longphase haplotagged.bam')

parser.add_argument('-f', '--bam', required=True, help="cellranger bam")
parser.add_argument('-o', '--out', required=True, help="output fastq name")
args = parser.parse_args()

fn = args.bam  # "haplotagged.bam"#files[0]
bam = pysam.AlignmentFile(fn, "rb")

with open(args.out, 'w') as fastq:
    for (index, read) in enumerate(bam):
        if read.has_tag('HP'):
            haplotype = read.get_tag('HP')
            phaseset = read.get_tag('PS')
        if read.is_secondary or read.is_supplementary or read.seq is None:
            continue

        readname = read.query_name
        if read.has_tag('HP'):
            fastq.write("@" + read.query_name + ";HP:" + str(haplotype) + ";PS:" + str(phaseset) + "\n")
        else:
            fastq.write("@" + read.query_name + "\n")
        fastq.write(read.seq + "\n")
        fastq.write("+\n")
        fastq.write(read.qual + "\n")
