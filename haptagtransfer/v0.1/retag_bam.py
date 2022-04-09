#!/usr/bin/env python

import pysam
import argparse


parser = argparse.ArgumentParser(description = 'Retag reads with their haplotype')

parser.add_argument('-s', '--sam', required = True, help = "sam file")
parser.add_argument('-o', '--out', required = True, help = "output (will be a bam)")
args = parser.parse_args()

bam = pysam.AlignmentFile(args.sam)

bamout = pysam.AlignmentFile(args.out,'wb', template = bam)

for read in bam:
    qname = read.qname
    tokens = qname.split(";")
    if ( len(tokens) == 3 ):
        read.set_tag("HP", int(tokens[-2][-1]))
        read.set_tag("PS", int(tokens[-1][3:]))
    bamout.write(read)

bamout.close()