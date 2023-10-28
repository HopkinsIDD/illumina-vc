#!/bin/bash
set -e 
set -o pipefail

#run align_reference.sh to obtain merged bam file before running this
#use samtools to create depth bedfile
#use bcftools to call variants and create consensus genome


# ---------------- Input Paramters and Files  ---------------- #

# sample name (to be run in parallel for multiple samples)
SAMPLENAME=$1

# path to project directory (should be consistent for all scripts)
DIR=$2
DIR=${DIR%/} # remove trailing slash from directory name if necessary

# path to input bam file
# If you have run 'align_reference.sh' 
# This file should be in 02_merged and with a suffix '.merged.mkdup.bam'
BAM=$3

# number of available cores
NUM_CORES=$4

# path to reference genome
REF_GENOME=$5


# set desired variant calling parameters
SNP_SUPPORT=$6 # required percent(0-1) of mapped reads to support alternate allele
REQ_DEPTH=$7 # required at least total depth
REQ_STRAND_DEPTH=$8 # required at least the number of reads matching alternate allele per strand


cd $DIR
echo $SAMPLENAME

# create folder to store tmp files for checking file completion

if [ ! -d check_tmp_variant ]; then
	mkdir check_tmp_variant
fi


# ---------------- Check bam file existence exists  ---------------- # 


if [ ! -f $BAM ]; then
	echo "ERROR: no merged bam file with marked duplicates"
	exit
fi

BAM_DIR=${BAM%/*}


# ---------------- Call variants and build consensus genome ---------------- #

date

# index reference with faidx
samtools faidx $REF_GENOME

echo "CALLING VARIANTS"

# ensure output directory exists
if [ ! -d $DIR/03_variants ]; then
	mkdir $DIR/03_variants
fi

if [ ! -d $DIR/04_assembly ]; then
	mkdir $DIR/04_assembly
fi

# create bed file from bam file
if [ ! -f $BAM_DIR/$SAMPLENAME.depth ];then
	samtools depth -a $BAM -q 20 -Q 30 > \
	$BAM_DIR/$SAMPLENAME.depth
	touch check_tmp_variant/.$SAMPLENAME.bam_depth
elif [ ! -e check_tmp_variant/.$SAMPLENAME.bam_depth ]; then
	rm $BAM_DIR/$SAMPLENAME.depth
	samtools depth -a $BAM -q 20 -Q 30 > \
	$BAM_DIR/$SAMPLENAME.depth
	touch check_tmp_variant/.$SAMPLENAME.bam_depth
fi


# filter bed file for sites with less than desired coverage

if [ ! -f $DIR/04_assembly/$SAMPLENAME.depthmask.bed ]; then
	awk -v depth="$REQ_DEPTH" '$3 < depth {printf "%s\t%d\t%d\n", $1, $2 - 1, $2}' \
	$BAM_DIR/$SAMPLENAME.depth > $DIR/04_assembly/$SAMPLENAME.depthmask.bed
	touch check_tmp_variant/.$SAMPLENAME.bam_depth_low
elif [ ! -e check_tmp_variant/.$SAMPLENAME.bam_depth_low ]; then
	rm $DIR/04_assembly/$SAMPLENAME.depthmask.bed
	awk -v depth="$REQ_DEPTH" '$3 < depth {printf "%s\t%d\t%d\n", $1, $2 - 1, $2}' \
	$BAM_DIR/$SAMPLENAME.depth > $DIR/04_assembly/$SAMPLENAME.depthmask.bed
	touch check_tmp_variant/.$SAMPLENAME.bam_depth_low
fi


# get mpileup of reads at each position and call variants with bcftools
# this will take several minutes
if [ ! -f $DIR/03_variants/$SAMPLENAME.bcftools.vcf ]; then
	bcftools mpileup --threads $NUM_CORES -d 2000 -q 30 -Q 20 -B -a INFO/AD,INFO/ADF,INFO/ADR -Ou \
	-f $REF_GENOME $BAM | \
	bcftools call --threads $NUM_CORES -mv -Ov --ploidy 1 -o $DIR/03_variants/$SAMPLENAME.bcftools.vcf
	touch check_tmp_variant/.$SAMPLENAME.vcf_bcftools
elif [ ! -e check_tmp_variant/.$SAMPLENAME.vcf_bcftools ]; then
	rm $DIR/03_variants/$SAMPLENAME.bcftools.vcf
	bcftools mpileup --threads $NUM_CORES -d 2000 -q 30 -Q 20 -B -a INFO/AD,INFO/ADF,INFO/ADR -Ou \
	-f $REF_GENOME $BAM | \
	bcftools call --threads $NUM_CORES -mv -Ov --ploidy 1 -o $DIR/03_variants/$SAMPLENAME.bcftools.vcf
	touch check_tmp_variant/.$SAMPLENAME.vcf_bcftools
fi

# filter variants on hard coded support and read depth

if [ ! -f $DIR/03_variants/$SAMPLENAME.bcftools.filt.vcf ]; then
	bcftools filter --no-version -i "INFO/AD[1]>$REQ_DEPTH" $DIR/03_variants/$SAMPLENAME.bcftools.vcf | \
	bcftools filter --no-version -i "(INFO/AD[1])/(INFO/AD[0]+INFO/AD[1])>$SNP_SUPPORT" | \
	bcftools filter --no-version -i "INFO/ADF[1]>$REQ_STRAND_DEPTH && INFO/ADR[1]>$REQ_STRAND_DEPTH" \
	-o $DIR/03_variants/$SAMPLENAME.bcftools.filt.vcf
	touch check_tmp_variant/.$SAMPLENAME.vcf_bcftools_filt
elif [ ! -e check_tmp_variant/.$SAMPLENAME.vcf_bcftools_filt ]; then
	rm $DIR/03_variants/$SAMPLENAME.bcftools.filt.vcf
	bcftools filter --no-version -i "INFO/AD[1]>$REQ_DEPTH" $DIR/03_variants/$SAMPLENAME.bcftools.vcf | \
	bcftools filter --no-version -i "(INFO/AD[1])/(INFO/AD[0]+INFO/AD[1])>$SNP_SUPPORT" | \
	bcftools filter --no-version -i "INFO/ADF[1]>$REQ_STRAND_DEPTH && INFO/ADR[1]>$REQ_STRAND_DEPTH" \
	-o $DIR/03_variants/$SAMPLENAME.bcftools.filt.vcf
	touch check_tmp_variant/.$SAMPLENAME.vcf_bcftools_filt
fi

# mask bcftool filtered positions
# concat bcftool filtered-masked positions with depth-masked bed file

if [ ! -f $DIR/04_assembly/$SAMPLENAME.allmask.bed ]; then
	bcftools filter --no-version -i "(INFO/AD[1])/(INFO/AD[0]+INFO/AD[1]) < $SNP_SUPPORT" \
	$DIR/03_variants/$SAMPLENAME.bcftools.vcf -o $DIR/03_variants/$SAMPLENAME.bcftools.lowsupport.vcf
	awk '(/^[^#]/ && length($4) == length($5)) {printf "%s\t%d\t%d\n", $1, $2 - 1, $2}' \
	$DIR/03_variants/$SAMPLENAME.bcftools.lowsupport.vcf > $DIR/03_variants/$SAMPLENAME.lowsupport.bed

	bcftools filter --no-version -i "INFO/ADF[1] < $REQ_STRAND_DEPTH || INFO/ADR[1] < $REQ_STRAND_DEPTH" \
	$DIR/03_variants/$SAMPLENAME.bcftools.vcf -o $DIR/03_variants/$SAMPLENAME.bcftools.lowstranddepth.vcf
	awk '(/^[^#]/ && length($4) == length($5)) {printf "%s\t%d\t%d\n", $1, $2 - 1, $2}' \
	$DIR/03_variants/$SAMPLENAME.bcftools.lowstranddepth.vcf > $DIR/03_variants/$SAMPLENAME.lowstranddepth.bed


	cat $DIR/03_variants/$SAMPLENAME.lowsupport.bed $DIR/03_variants/$SAMPLENAME.lowstranddepth.bed > \
	$DIR/03_variants/$SAMPLENAME.vcffiltered.bed

	cat $DIR/03_variants/$SAMPLENAME.vcffiltered.bed $DIR/04_assembly/$SAMPLENAME.depthmask.bed | \
	uniq > $DIR/04_assembly/$SAMPLENAME.allmask.bed

	touch check_tmp_variant/.$SAMPLENAME.all_mask_bed

elif [ ! -e check_tmp_variant/.$SAMPLENAME.all_mask_bed ]; then

	rm -f $DIR/03_variants/$SAMPLENAME.bcftools.lowsupport.vcf
	rm -f $DIR/03_variants/$SAMPLENAME.lowsupport.bed
	rm -f $DIR/03_variants/$SAMPLENAME.bcftools.lowstranddepth.vcf
	rm -f $DIR/03_variants/$SAMPLENAME.lowstranddepth.bam_depth
	rm -f $DIR/03_variants/$SAMPLENAME.vcffiltered.bed
	rm -f $DIR/04_assembly/$SAMPLENAME.allmask.bed

	bcftools filter --no-version -i "(INFO/AD[1])/(INFO/AD[0]+INFO/AD[1]) < $SNP_SUPPORT" \
	$DIR/03_variants/$SAMPLENAME.bcftools.vcf -o $DIR/03_variants/$SAMPLENAME.bcftools.lowsupport.vcf
	awk '(/^[^#]/ && length($4) == length($5)) {printf "%s\t%d\t%d\n", $1, $2 - 1, $2}' \
	$DIR/03_variants/$SAMPLENAME.bcftools.lowsupport.vcf > $DIR/03_variants/$SAMPLENAME.lowsupport.bed

	bcftools filter --no-version -i "INFO/ADF[1] < $REQ_STRAND_DEPTH || INFO/ADR[1] < $REQ_STRAND_DEPTH" \
	$DIR/03_variants/$SAMPLENAME.bcftools.vcf -o $DIR/03_variants/$SAMPLENAME.bcftools.lowstranddepth.vcf
	awk '(/^[^#]/ && length($4) == length($5)) {printf "%s\t%d\t%d\n", $1, $2 - 1, $2}' \
	$DIR/03_variants/$SAMPLENAME.bcftools.lowstranddepth.vcf > $DIR/03_variants/$SAMPLENAME.lowstranddepth.bed


	cat $DIR/03_variants/$SAMPLENAME.lowsupport.bed $DIR/03_variants/$SAMPLENAME.lowstranddepth.bed > \
	$DIR/03_variants/$SAMPLENAME.vcffiltered.bed

	cat $DIR/03_variants/$SAMPLENAME.vcffiltered.bed $DIR/04_assembly/$SAMPLENAME.depthmask.bed | \
	uniq > $DIR/04_assembly/$SAMPLENAME.allmask.bed

	touch check_tmp_variant/.$SAMPLENAME.all_mask_bed
fi

# left align and normalize indels
# then remove insertions to avoid issues with multi-alignment later on

if [ ! -f $DIR/03_variants/$SAMPLENAME.bcftools.filt.norm.vcf.gz ]; then
	bcftools norm --no-version -f $REF_GENOME $DIR/03_variants/$SAMPLENAME.bcftools.filt.vcf | \
	bcftools filter --no-version --exclude 'strlen(REF)<strlen(ALT)' -Oz \
	-o $DIR/03_variants/$SAMPLENAME.bcftools.filt.norm.vcf.gz
	touch check_tmp_variant/.$SAMPLENAME.vcf_bcftools_filt_norm
elif [ ! -e check_tmp_variant/.$SAMPLENAME.vcf_bcftools_filt_norm ]; then
	rm $DIR/03_variants/$SAMPLENAME.bcftools.filt.norm.vcf.gz
	bcftools norm --no-version -f $REF_GENOME $DIR/03_variants/$SAMPLENAME.bcftools.filt.vcf | \
	bcftools filter --no-version --exclude 'strlen(REF)<strlen(ALT)' -Oz \
	-o $DIR/03_variants/$SAMPLENAME.bcftools.filt.norm.vcf.gz
	touch check_tmp_variant/.$SAMPLENAME.vcf_bcftools_filt_norm
fi	

echo "CREATING CONSENSUS GENOME"

# apply variants to create consensus sequence
# mask sites with less than desired coverage

if [ ! -f $DIR/04_assembly/$SAMPLENAME.bcftools.consensus.fasta ]; then
	bcftools index --threads $NUM_CORES $DIR/03_variants/$SAMPLENAME.bcftools.filt.norm.vcf.gz
	cat $REF_GENOME | bcftools consensus --mark-del N -m $DIR/04_assembly/$SAMPLENAME.depthmask.bed \
	$DIR/03_variants/$SAMPLENAME.bcftools.filt.norm.vcf.gz > \
	$DIR/04_assembly/$SAMPLENAME.bcftools.consensus.fasta
	touch check_tmp_variant/.$SAMPLENAME.bcftools_consensus_fasta
elif [ ! -e check_tmp_variant/.$SAMPLENAME.bcftools_consensus_fasta ]; then
	rm $DIR/04_assembly/$SAMPLENAME.bcftools.consensus.fasta
	bcftools index --threads $NUM_CORES $DIR/03_variants/$SAMPLENAME.bcftools.filt.norm.vcf.gz
	cat $REF_GENOME | bcftools consensus --mark-del N -m $DIR/04_assembly/$SAMPLENAME.depthmask.bed \
	$DIR/03_variants/$SAMPLENAME.bcftools.filt.norm.vcf.gz > \
	$DIR/04_assembly/$SAMPLENAME.bcftools.consensus.fasta
	touch check_tmp_variant/.$SAMPLENAME.bcftools_consensus
fi


# merge two chromosomes to create single fasta
# update the fasta header to indicate the correct sample

header=">$SAMPLENAME"

if [ ! -f $DIR/04_assembly/$SAMPLENAME.bcftools.fasta ]; then
	union -filter $DIR/04_assembly/$SAMPLENAME.bcftools.consensus.fasta > $DIR/04_assembly/$SAMPLENAME.bcftools.fasta
	sed -i "1s/.*/$header/" $DIR/04_assembly/$SAMPLENAME.bcftools.fasta
	touch check_tmp_variant/.$SAMPLENAME.bcftools_merge_consensus
elif [ ! -e check_tmp_variant/.$SAMPLENAME.bcftools_merge_consensus ]; then
	rm $DIR/04_assembly/$SAMPLENAME.bcftools.fasta
	union -filter $DIR/04_assembly/$SAMPLENAME.bcftools.consensus.fasta > $DIR/04_assembly/$SAMPLENAME.bcftools.fasta
	sed -i "1s/.*/$header/" $DIR/04_assembly/$SAMPLENAME.bcftools.fasta
	touch check_tmp_variant/.$SAMPLENAME.bcftools_merge_consensus
fi



# check length of sequences
RESULT_LENGTH=$(awk '/^>/ {if (seqlen){print seqlen};seqlen=0;next; } \
				{ seqlen += length($0)}END{print seqlen}' \
				$DIR/04_assembly/$SAMPLENAME.bcftools.fasta)

if [ $RESULT_LENGTH -eq 4033501 ]; then
	echo "reference-based consensus created successfully"
else
	echo "WARNING: sequence length didn't equal to reference genome"
fi


date

