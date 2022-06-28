#!/bin/bash
set -e 
set -o pipefail

#Use bwa to index reference genome
#then use bwa, samtools to align reads to reference genome
#use samtools to sort and merge bam files if more than one lane
#use picard to mark duplicate for merged bam file

# ---------------- Input Paramters and Files  ---------------- #

# path to project directory (should be consistent for all scripts)
DIR=$1
DIR=${DIR%/} # remove trailing slash from directory name if necessary

# path to FASTQ data
FASTQ=$2
FASTQ=${FASTQ%/} # remove trailing slash from directory name if necessary

# sample name 
SAMPLENAME=$3

# indicate number of lanes to merge data over
LANES=$4

# number of available cores
NUM_CORES=$5

# path to picard jar file
PICARD_PATH=$6 

# path to reference genome
REF_GENOME=$7


cd $DIR
echo $SAMPLENAME



# create folder to store tmp files for checking file completion

if [ ! -d check_tmp_align ]; then
	mkdir check_tmp_align
fi

# ---------------- Map reads to reference ---------------- #

date

# index reference genome
if [ ! -f $REF_GENOME.ann ]; then
	
	bwa index $REF_GENOME
	touch check_tmp_align/.index_ref_check

elif [ ! -e check_tmp_align/.index_ref_check ]; then

	rm $REF_GENOME.ann
	bwa index $REF_GENOME
	touch check_tmp_align/.index_ref_check

fi

echo "MAPPING READS TO REFERENCE"

# set up array of all lanes to loop through
declare -a LARRAY
for i in $(seq 1 $LANES); do LARRAY+=("L00"$i); done

# loop through all lanes
for i in "${LARRAY[@]}"; do

	LANE=$i

	# ensure output directory exists
	if [ ! -d $DIR/01_mapped ]; then
		mkdir $DIR/01_mapped
	fi

	# align paired end reads to reference genome (required indexed reference)
	# this will take a few minutes
	if [ ! -f $DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.sam ]; then
		bwa mem -M -t $NUM_CORES $REF_GENOME $FASTQ/$SAMPLENAME\_$LANE\_R1_001.fastq.gz \
		$FASTQ/$SAMPLENAME\_$LANE\_R2_001.fastq.gz > $DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.sam
		touch check_tmp_align/.$SAMPLENAME\_$LANE_map_sam_check
	elif [ ! -e check_tmp_align/.$SAMPLENAME\_$LANE_map_sam_check ]; then
		rm $DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.sam
		bwa mem -M -t $NUM_CORES $REF_GENOME $FASTQ/$SAMPLENAME\_$LANE\_R1_001.fastq.gz \
		$FASTQ/$SAMPLENAME\_$LANE\_R2_001.fastq.gz > $DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.sam
		touch check_tmp_align/.$SAMPLENAME\_$LANE_map_sam_check
	fi

	# convert SAM to BAM and fix read pairing information and flags
	# this will take a minute
	if [ ! -f $DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.bam ]; then
		samtools fixmate -O bam $DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.sam \
		$DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.bam
		touch check_tmp_align/.$SAMPLENAME\_$LANE_sam_bam
	elif [ ! -e check_tmp_align/.$SAMPLENAME\_$LANE_sam_bam ]; then
		rm $DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.bam
		samtools fixmate -O bam $DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.sam \
		$DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.bam
		touch check_tmp_align/.$SAMPLENAME\_$LANE_sam_bam
	fi


	# sort BAM file
	if [ ! -f $DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.sorted.bam ]; then
		samtools sort -T $DIR/01_mapped/$SAMPLENAME\_$LANE.mapped -O bam -@ $NUM_CORES \
		-o $DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.sorted.bam \
		$DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.bam
		touch check_tmp_align/.$SAMPLENAME\_$LANE_sort_bam
	elif [ ! -e check_tmp_align/.$SAMPLENAME\_$LANE_sort_bam ]; then
		rm $DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.sorted.bam
		samtools sort -T $DIR/01_mapped/$SAMPLENAME\_$LANE.mapped -O bam -@ $NUM_CORES \
		-o $DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.sorted.bam \
		$DIR/01_mapped/$SAMPLENAME\_$LANE.mapped.bam
		touch check_tmp_align/.$SAMPLENAME\_$LANE_sort_bam
	fi

done

# create list of files to merge from different lanes
if [ ! -f $DIR/01_mapped/$SAMPLENAME.lanebams.txt ]; then
	ls $DIR/01_mapped/$SAMPLENAME\_L00*.mapped.sorted.bam >> \
	$DIR/01_mapped/$SAMPLENAME.lanebams.txt
else
	rm $DIR/01_mapped/$SAMPLENAME.lanebams.txt
	ls $DIR/01_mapped/$SAMPLENAME\_L00*.mapped.sorted.bam >> \
	$DIR/01_mapped/$SAMPLENAME.lanebams.txt
fi

# ensure output directory exists
if [ ! -d $DIR/02_merged ]; then
	mkdir $DIR/02_merged
fi

# merge data from multiple lanes if needed
# this will take a minute
if [ "$LANES" -gt 1 ]; then
	if [ ! -f $DIR/02_merged/$SAMPLENAME.merged.bam ]; then
		samtools merge -f $DIR/02_merged/$SAMPLENAME.merged.bam \
		-b $DIR/01_mapped/$SAMPLENAME.lanebams.txt
		touch check_tmp_align/.$SAMPLENAME\_merge_bam
	elif [ ! -e check_tmp_align/.$SAMPLENAME\_merge_bam ]; then
		rm $DIR/02_merged/$SAMPLENAME.merged.bam
		samtools merge -f $DIR/02_merged/$SAMPLENAME.merged.bam \
		-b $DIR/01_mapped/$SAMPLENAME.lanebams.txt
		touch check_tmp_align/.$SAMPLENAME\_merge_bam
	fi
else
	if [ ! -f $DIR/02_merged/$SAMPLENAME.merged.bam ]; then
		cp $DIR/01_mapped/$SAMPLENAME\_L001.mapped.sorted.bam \
		$DIR/02_merged/$SAMPLENAME.merged.bam
		touch check_tmp_align/.$SAMPLENAME\_merge_bam
	elif [ ! -e check_tmp_align/.$SAMPLENAME\_merge_bam ]; then
		rm $DIR/02_merged/$SAMPLENAME.merged.bam
		cp $DIR/01_mapped/$SAMPLENAME\_L001.mapped.sorted.bam \
		$DIR/02_merged/$SAMPLENAME.merged.bam
		touch check_tmp_align/.$SAMPLENAME\_merge_bam
	fi
fi

# mark duplicates in merged bam
if [ ! -f $DIR/02_merged/$SAMPLENAME.merged.mkdup.bam ]; then
	java -Xmx16g -jar $PICARD_PATH MarkDuplicates ASSUME_SORTED=true \
	INPUT=$DIR/02_merged/$SAMPLENAME.merged.bam OUTPUT=$DIR/02_merged/$SAMPLENAME.merged.mkdup.bam \
	METRICS_FILE=$DIR/02_merged/$SAMPLENAME.mkdup.metrics.out
	touch check_tmp_align/.$SAMPLENAME\_merge_mkdup_bam
elif [ ! -e check_tmp_align/.$SAMPLENAME\_merge_mkdup_bam ]; then
	rm $DIR/02_merged/$SAMPLENAME.merged.mkdup.bam
	rm $DIR/02_merged/$SAMPLENAME.mkdup.metrics.out
	java -Xmx16g -jar $PICARD_PATH MarkDuplicates ASSUME_SORTED=true \
	INPUT=$DIR/02_merged/$SAMPLENAME.merged.bam OUTPUT=$DIR/02_merged/$SAMPLENAME.merged.mkdup.bam \
	METRICS_FILE=$DIR/02_merged/$SAMPLENAME.mkdup.metrics.out
	touch check_tmp_align/.$SAMPLENAME\_merge_mkdup_bam
fi


# calculate alignment statistics for merged bam (optional)
# if [ ! -f $DIR/02_merged/$SAMPLENAME.merged.mkdup.stats.out ]; then
# 	samtools flagstat $DIR/02_merged/$SAMPLENAME.merged.mkdup.bam > \
# 	$DIR/02_merged/$SAMPLENAME.merged.mkdup.stats.out
# 	touch check_tmp_align/.$SAMPLENAME\_merge_mkdup_stat
# elif [ ! -e check_tmp_align/.$SAMPLENAME\_merge_mkdup_stat ]; then
# 	rm $DIR/02_merged/$SAMPLENAME.merged.mkdup.stats.out
# 	samtools flagstat $DIR/02_merged/$SAMPLENAME.merged.mkdup.bam > \
# 	$DIR/02_merged/$SAMPLENAME.merged.mkdup.stats.out
# 	touch check_tmp_align/.$SAMPLENAME\_merge_mkdup_stat
# fi

printf "%s mapped to refernce genome" "$SAMPLENAME"
date