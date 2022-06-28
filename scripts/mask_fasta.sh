#! /bin/bash
set -e 
set -o pipefail

# Mask sequences with a given masking file
#use bedtools to mask sites with known recombination

# ---------------- Input Paramters and Files  ---------------- #

# path to project directory (should be consistent for all scripts)
DIR=$1
DIR=${DIR%/} # remove trailing slash from directory name if necessary

# sample name (to be run in parallel for multiple samples)
SAMPLENAME=$2

# path to assembled genome fasta file
# If you have run 'align_reference.sh' and 'call_variant.sh'
# This file should be in 04_assembly and with a suffix '.bcftools.fasta'
FASTA=$3

# path to gff file with recombinant sites to mask
RECOMBMASK=$4

# chrom name in provided gff file
chrom=$5

cd $DIR
echo $SAMPLENAME

# create folder to store tmp files for checking file completion

if [ ! -d check_tmp_mask ]; then
	mkdir check_tmp_mask
fi

echo "Mask sites with the given masking file"

# ---------------- Mask sites with the given masking file  ---------------- #

date

header=$(head -1 $FASTA)
FASTA_DIR=${FASTA%/*}

if [ ! -f $FASTA_DIR/$SAMPLENAME.mask.fasta ]; then
	sed -i "1s/.*/$chrom/" $FASTA # temporarily change fasta header to match gff
	bedtools maskfasta -fi $FASTA \
	-bed $RECOMBMASK -fo $FASTA_DIR/$SAMPLENAME.mask.fasta
	sed -i "1s/.*/$header/" $FASTA # change fasta header back to sample name
	sed -i "1s/.*/$header/" $FASTA_DIR/$SAMPLENAME.mask.fasta # update header in new file
	touch check_tmp_mask/.$SAMPLENAME.recomb_mask_consensus
elif [ ! -e check_tmp_mask/.$SAMPLENAME.recomb_mask_consensus ]; then
	rm $FASTA_DIR/$SAMPLENAME.mask.fasta
	sed -i "1s/.*/$chrom/" $FASTA # temporarily change fasta header to match gff
	bedtools maskfasta -fi $FASTA \
	-bed $RECOMBMASK -fo $FASTA_DIR/$SAMPLENAME.mask.fasta
	sed -i "1s/.*/$header/" $FASTA # change fasta header back to sample name
	sed -i "1s/.*/$header/" $FASTA_DIR/$SAMPLENAME.mask.fasta # update header in new file
	touch check_tmp_mask/.$SAMPLENAME.recomb_mask_consensus
fi	

echo "Masking completed"
date