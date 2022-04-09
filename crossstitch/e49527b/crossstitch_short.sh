#!/bin/bash

#set -xv
set -e

if [ $# -ne 5 ]
then
  echo "USAGE: crossstitch.sh phased_snps.vcf unphased_structural_variants.vcf long_reads.bam genome.fa outputprefix karyotype refine"
  echo ""
  echo "Details:"
  echo "  phased_snps.vcf:                   VCF file of phased SNP and indel variants. Recommend LongRanger (10X only) or HapCUT2 (HiC and/or 10X)"
  echo "  phased_structural_variants.vcf:    VCF file of structural variants identified using Sniffles"
  echo "  genome.fa:                         Reference genome used"
  echo "  outputprefix:                      Prefix for output files"
  echo "  karyotype:                         "xy" or "xx", used to ensure sex chromosomes are correctly used"
  exit
fi

BINDIR=`cd "$(dirname "$0")" ; pwd -P`
PHASEDSNPS=$1
STRUCTURALVARIANTS=$2
GENOME=$3
OUTPREFIX=$4
KARYOTYPE=$5

VCF2DIPLOIDJAR=$BINDIR/../vcf2diploid/vcf2diploid.jar
EXTRACTHAIRS=extractHAIRS
GZIP=pigz

echo "crossstich.sh"
echo "  BINDIR: $BINDIR"
echo "  PHASEDSNPS: $PHASEDSNPS"
echo "  STRUCTURALVARIANTS: $STRUCTURALVARIANTS"
echo "  GENOME: $GENOME"
echo "  KARYOTYPE: $KARYOTYPE"
echo 
echo "  OUT: $OUTPREFIX"
echo
echo

## Sanity check parameters

if [[ $KARYOTYPE != "xy" && $KARYOTYPE != "xx" ]]
then
  echo "Unknown karyotype: $KARYOTYPE (must be xy or xx)"
  exit
fi

if [ ! -r $GENOME ]
then
  echo "Cannot read genome file: $GENOME"
  exit
fi

if [ ! -r $STRUCTURALVARIANTS ]
then
  echo "Cannot read sv vcf file: $STRUCTURALVARIANTS"
  exit
fi

if [ ! -r $PHASEDSNPS ]
then
  echo "Cannot read phased snp file: $PHASEDSNPS"
  exit
fi

## Sanity checks passed, begin analysis

GENOME="$(cd "$(dirname "$GENOME")" ; pwd -P)/$(basename $GENOME)"
AS=$OUTPREFIX.alleleseq
VCFID=`head -5000 $PHASEDSNPS | grep '#CHROM' | awk '{print $10}'`

if [ ! -r $OUTPREFIX.scrubbed.vcf ]
then
  echo "Scrubbing SV calls"
  (java -cp $BINDIR RemoveInvalidVariants $STRUCTURALVARIANTS $OUTPREFIX.scrubbed.vcf) >& $OUTPREFIX.scrubbed.log
fi

echo $OUTPREFIX > samples
bcftools reheader -s samples $PHASEDSNPS | bcftools view -o $OUTPREFIX.scrubbed.reheader.snv.indel.vcf.gz 
bcftools reheader -s samples $OUTPREFIX.scrubbed.vcf | bcftools view -o $OUTPREFIX.scrubbed.reheader.svs.vcf.gz 
bcftools index $OUTPREFIX.scrubbed.reheader.snv.indel.vcf.gz 
bcftools index $OUTPREFIX.scrubbed.reheader.svs.vcf.gz
bcftools concat -a -o $OUTPREFIX.spliced.vcf $OUTPREFIX.scrubbed.reheader.snv.indel.vcf.gz $OUTPREFIX.scrubbed.reheader.svs.vcf.gz

if [ ! -r $OUTPREFIX.spliced.scrubbed.vcf ]
then
  echo "Final scrub to remove overlapping spliced variants"
  (java -cp $BINDIR RemoveInvalidVariants -o $KARYOTYPE $OUTPREFIX.spliced.vcf $OUTPREFIX.spliced.scrubbed.vcf) >& $OUTPREFIX.spliced.scrubbed.log
fi

if [ ! -r $OUTPREFIX.spliced.scrubbed.vcf.gz ]
then
  echo "compressing spliced scrubbed vcf"
  $GZIP -c $OUTPREFIX.spliced.scrubbed.vcf > $OUTPREFIX.spliced.scrubbed.vcf.gz
fi

if [ ! -r $AS ]
then
  mkdir -p $AS
  pushd $AS

  echo "constructing diploid sequence with SNPs and SVs"
  java -Xmx400000m -jar $VCF2DIPLOIDJAR -id $VCFID -pass -chr $GENOME -vcf ../$OUTPREFIX.spliced.scrubbed.vcf >& vcf2diploid.log
  popd
fi

# if [ ! -r $AS.raw.tgz ]
# then
#   echo "tarring up alleleseq"
#   tar czvf $AS.raw.tgz $AS
# fi

if [ ! -r $AS/raw/ ]
then
  echo "renaming files"

  pushd $AS

  mkdir -p raw
  mkdir -p raw/attic

  for i in `ls *_$VCFID*`
  do 
    echo " $i"
    j=`echo $i | sed "s/_$VCFID//"`
    mv $i $OUTPREFIX.$j
  done

  for i in `/bin/ls *_paternal*`
  do 
    j=`echo $i | sed s/_paternal/.hap1/`
    mv $i $j
  done

  for i in `/bin/ls *_maternal*`
  do 
    j=`echo $i | sed s/_maternal/.hap2/`
    mv $i $j
  done

  mv paternal.chain hap1.chain
  mv maternal.chain hap2.chain

  mv hap1.chain raw/$OUTPREFIX.hap1.chain
  mv hap2.chain raw/$OUTPREFIX.hap2.chain

  if [ -r $OUTPREFIX.chrM.hap2.fa ]
  then
    mv $OUTPREFIX.chrM.hap2.fa raw/attic
  fi

  if [[ $KARYOTYPE == "xy" ]]
  then
    echo "xy sample, making X and Y haploid"

    if [ -r $OUTPREFIX.chrX.hap2.fa ]
    then
      mv *chrX.hap2.fa *chrY.hap2.fa raw/attic
    fi

    cat raw/$OUTPREFIX.hap1.chain                                   | sed 's/paternal/hap1/' > $OUTPREFIX.hap1.chain
    $BINDIR/removechain.pl raw/$OUTPREFIX.hap2.chain chrM chrX chrY | sed 's/maternal/hap2/' > $OUTPREFIX.hap2.chain
  else
    echo "xx sample, stashing Y chromosome"

    if [ -r $OUTPREFIX.chrY.hap1.fa ]
    then
      mv *chrY* raw/attic
    fi

    $BINDIR/removechain.pl raw/$OUTPREFIX.hap1.chain chrY           | sed 's/paternal/hap1/' > $OUTPREFIX.hap1.chain
    $BINDIR/removechain.pl raw/$OUTPREFIX.hap2.chain chrM chrY      | sed 's/maternal/hap2/' > $OUTPREFIX.hap2.chain
  fi

  echo "fixing map files"
  for i in `/bin/ls *.map`
  do
    echo " $i"
    mv $i raw/
    sed 's/PAT/HAP1/' raw/$i | sed 's/MAT/HAP2/' > $i
  done

  echo "fixing chromosome files"
  for i in `/bin/ls *hap1.fa`
  do
    echo " $i"
    mv $i raw/
    sed 's/paternal/hap1/' raw/$i > $i
  done

  for i in `/bin/ls *hap2.fa`
  do
    echo " $i"
    mv $i raw/
    sed 's/maternal/hap2/' raw/$i > $i
  done

  popd
fi

if [ ! -r $AS/$OUTPREFIX.hap1.fa.gz ]
then
  echo "Making final diploid genome"
  cat `ls $AS/*.chr*.hap1.fa | sort -V` > $AS/$OUTPREFIX.hap1.fa &
  cat `ls $AS/*.chr*.hap2.fa | sort -V` > $AS/$OUTPREFIX.hap2.fa &

  wait

  echo "compressing genome files"
  $GZIP -c $AS/$OUTPREFIX.hap1.fa > $AS/$OUTPREFIX.hap1.fa.gz
  $GZIP -c $AS/$OUTPREFIX.hap2.fa > $AS/$OUTPREFIX.hap2.fa.gz
fi