#!/usr/bin/env python

import pysam
from Bio.Seq import Seq
import argparse


parser = argparse.ArgumentParser(description = 'Add ALT allele sequence for called Inversions')

parser.add_argument('-v', '--vcf', required = True, help = "Input vcf")
parser.add_argument('-f', '--fasta', required = True, help = "reference genome fasta")
parser.add_argument('-o', '--out', required = True, help = "output vcf")
args = parser.parse_args()

vcf_in = pysam.VariantFile(args.vcf)
vcf_out = pysam.VariantFile(args.out, 'w', header=vcf_in.header)

refseq = pysam.FastaFile(args.fasta)

for rec in vcf_in.fetch():
    if ( rec.alts[0] == "<INV>" ):
        # print(rec.alts)
        # print(rec.info['SVLEN'])
        ref = refseq.fetch(rec.chrom, rec.pos - 1, rec.pos + rec.info['SVLEN'] - 1)
        # print(len(str(ref)))
        rec.alts = (str(Seq(ref).reverse_complement()), )
        rec.ref = "X" * len(ref)
        # print(rec.alts)
    vcf_out.write(rec)

