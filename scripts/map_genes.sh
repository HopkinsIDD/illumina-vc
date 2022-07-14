#!/bin/bash
set -e 
set -o pipefail

# Align reads for a given sample to the following genes:
## ctxA -> cholera toxin; associated with pandemic lineages O1 and O139
## wbeO1 -> O1 antigen
## wbfO139 -> O139 antigen
## toxR -> virulence regulator gene
## tcpA classical
## tcpA El Tor

# Calculate metrics for each gene on mapped reads
## number and fraction of reads mapped
## depth of coverage
## percent coverage (positions with a minimum coverage depth)


# ---------------- Input Paramters and Files  ---------------- #

# sample name 
SAMPLENAME=$1

# path to project directory (create a new directory for typing purposes)
DIR=$2
DIR=${DIR%/} # remove trailing slash from directory name if necessary

# path to FASTQ data
FASTQ=$3
FASTQ=${FASTQ%/} # remove trailing slash from directory name if necessary

# indicate number of lanes to merge data over
LANES=$4

# number of available cores
NUM_CORES=$5

# path to typing reference genomes directory
REFDIR=$6
REFDIR=${REFDIR%/}

# path for output typing summary tables
OUTPUT=$7
OUTPUT=${OUTPUT%/}

cd $DIR
echo $SAMPLENAME


# set up array of all lanes to loop through
declare -a LARRAY
for i in $(seq 1 $LANES); do LARRAY+=("L00"$i); done


# ensure output directorys exists
if [ ! -d $DIR/typing_mapped ]; then
	mkdir $DIR/typing_mapped
fi

if [ ! -d $DIR/typing_merged ]; then
	mkdir $DIR/typing_merged
fi

if [ ! -d $OUTPUT ]; then
	mkdir $OUTPUT
fi

if [ ! -d $DIR/check_tmp_typing ]; then
	mkdir $DIR/check_tmp_typing
fi

#make sure the summary doesn't exist

if [ -f $OUTPUT/$SAMPLENAME\_typing_summary.txt ]; then
	rm $OUTPUT/$SAMPLENAME\_typing_summary.txt
fi


printf "SAMPLENAME\tgene\treads_mapped\treads\tfrac_mapped\tmean_depth\tpos_covered\ttotal_pos\tfrac_covered\n" >> \
	$OUTPUT/$SAMPLENAME\_typing_summary.txt


# ---------------- Map reads to reference and calculate metrics ---------------- #

# loop through all typing genes

date

for ref in $REFDIR/*.fasta; do

	# obtain the name of the gene
	gene=${ref%%.fasta}
	gene=${gene##*/}

	if [ ! -f $ref.ann ]; then
		bwa index $ref
		touch check_tmp_typing/.$gene\_ann_check
	elif [ ! -e check_tmp_typing/.$gene\_ann_check ]; then
		rm $ref.ann
		bwa index $ref
		touch check_tmp_typing/.$gene\_ann_check
	fi

	

	printf "MAPPING READS To %s" "$gene"

	# loop through lanes and align reads to reference
	for i in "${LARRAY[@]}"; do

		LANE=$i

		# align paired end reads to reference genome (required indexed reference)
		# this will take a few minutes
		if [ ! -f $DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.sam ]; then
			bwa mem -M -t $NUM_CORES $ref $FASTQ/$SAMPLENAME\_$LANE\_R1_001.fastq.gz \
			$FASTQ/$SAMPLENAME\_$LANE\_R2_001.fastq.gz > \
			$DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.sam
			touch check_tmp_typing/.$SAMPLENAME.$gene\_$LANE.sam.check
		elif [ ! -e check_tmp_typing/.$SAMPLENAME.$gene\_$LANE.sam.check ]; then
			rm $DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.sam
			bwa mem -M -t $NUM_CORES $ref $FASTQ/$SAMPLENAME\_$LANE\_R1_001.fastq.gz \
			$FASTQ/$SAMPLENAME\_$LANE\_R2_001.fastq.gz > \
			$DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.sam
			touch check_tmp_typing/.$SAMPLENAME.$gene\_$LANE.sam.check
		fi

		# convert SAM to BAM and fix read pairing information and flags
		# this will take a minute
		if [ ! -f $DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.bam ]; then
			samtools fixmate -O bam $DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.sam \
			$DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.bam
			touch check_tmp_typing/.$SAMPLENAME.$gene\_$LANE.bam.check
		elif [ ! -e check_tmp_typing/.$SAMPLENAME.$gene\_$LANE.bam.check ]; then
			rm $DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.bam
			samtools fixmate -O bam $DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.sam \
			$DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.bam
			touch check_tmp_typing/.$SAMPLENAME.$gene\_$LANE.bam.check			
		fi

		# sort BAM file
		if [ ! -f $DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.sorted.bam ]; then
			samtools sort -T $DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped -O bam \
			-@ $NUM_CORES -o $DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.sorted.bam \
			$DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.bam
			touch check_tmp_typing/.$SAMPLENAME.$gene\_$LANE.bam.sort.check
		elif [ ! -e check_tmp_typing/.$SAMPLENAME.$gene\_$LANE.bam.sort.check ]; then
			rm $DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.sorted.bam
			samtools sort -T $DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped -O bam \
			-@ $NUM_CORES -o $DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.sorted.bam \
			$DIR/typing_mapped/$SAMPLENAME.$gene\_$LANE.mapped.bam
			touch check_tmp_typing/.$SAMPLENAME.$gene\_$LANE.bam.sort.check
		fi

	done

	# create list of files to merge from different lanes
	if [ ! -f $DIR/typing_mapped/$SAMPLENAME.$gene.lanebams.txt ]; then
		ls $DIR/typing_mapped/$SAMPLENAME.$gene\_L00*.mapped.sorted.bam >> \
		$DIR/typing_mapped/$SAMPLENAME.$gene.lanebams.txt
	else
		rm $DIR/typing_mapped/$SAMPLENAME.$gene.lanebams.txt
		ls $DIR/typing_mapped/$SAMPLENAME.$gene\_L00*.mapped.sorted.bam >> \
		$DIR/typing_mapped/$SAMPLENAME.$gene.lanebams.txt
	fi



	# merge data from multiple lanes if needed
	# this will take a minute
	if [ "$LANES" -gt 1 ]; then
		if [ ! -f $DIR/typing_merged/$SAMPLENAME.$gene.merged.bam ]; then
			samtools merge -f $DIR/typing_merged/$SAMPLENAME.$gene.merged.bam \
			-b $DIR/typing_mapped/$SAMPLENAME.$gene.lanebams.txt
			touch check_tmp_typing/.$SAMPLENAME.$gene.bam.merge.check
		elif [ ! -e check_tmp_typing/.$SAMPLENAME.$gene.bam.merge.check ]; then
			rm $DIR/typing_merged/$SAMPLENAME.$gene.merged.bam
			samtools merge -f $DIR/typing_merged/$SAMPLENAME.$gene.merged.bam \
			-b $DIR/typing_mapped/$SAMPLENAME.$gene.lanebams.txt
			touch check_tmp_typing/.$SAMPLENAME.$gene.bam.merge.check
		fi
	else
		if [ ! -f $DIR/typing_merged/$SAMPLENAME.$gene.merged.bam ]; then
			cp $DIR/typing_mapped/$SAMPLENAME.$gene\_L001.mapped.sorted.bam \
			$DIR/typing_merged/$SAMPLENAME.$gene.merged.bam
			touch check_tmp_typing/.$SAMPLENAME.$gene.bam.merge.check
		elif [ ! -e check_tmp_typing/.$SAMPLENAME.$gene.bam.merge.check ]; then
			rm $DIR/typing_merged/$SAMPLENAME.$gene.merged.bam
			cp $DIR/typing_mapped/$SAMPLENAME.$gene\_L001.mapped.sorted.bam \
			$DIR/typing_merged/$SAMPLENAME.$gene.merged.bam
			touch check_tmp_typing/.$SAMPLENAME.$gene.bam.merge.check
		fi
	fi


	# index reference with faidx
	samtools faidx $ref

	# generate a depth file//
	if [ ! -f $DIR/typing_merged/$SAMPLENAME.$gene.merged.depth ]; then
		samtools depth -a $DIR/typing_merged/$SAMPLENAME.$gene.merged.bam > \
		$DIR/typing_merged/$SAMPLENAME.$gene.merged.depth
		touch check_tmp_typing/.$SAMPLENAME.$gene.depth.check
	elif [ ! -e check_tmp_typing/.$SAMPLENAME.$gene.depth.check ]; then
		rm $DIR/typing_merged/$SAMPLENAME.$gene.merged.depth
		samtools depth -a $DIR/typing_merged/$SAMPLENAME.$gene.merged.bam > \
		$DIR/typing_merged/$SAMPLENAME.$gene.merged.depth
		touch check_tmp_typing/.$SAMPLENAME.$gene.depth.check
	fi
	# calculate number of reads mapped
	reads_mapped=$(samtools view -c -F 4 $DIR/typing_merged/$SAMPLENAME.$gene.merged.bam)

	# calculate fraction of reads mapped
	# first calculate number of reads
	reads=$(samtools view -c $DIR/typing_merged/$SAMPLENAME.$gene.merged.bam)
	frac_mapped=$(echo "scale=7; $reads_mapped / $reads" | bc)

	# calculate mean depth
	mean_depth=$(awk '{ total += $3 } END { print total/NR }' \
		$DIR/typing_merged/$SAMPLENAME.$gene.merged.depth)

	# calculate percent of gene covered
	pos_covered=$(awk -v MIN_COV=$MIN_COV 'BEGIN {count = 0} $3 > MIN_COV {count++} END {print count}' \
		$DIR/typing_merged/$SAMPLENAME.$gene.merged.depth)
	total_pos=$(wc -l < $DIR/typing_merged/$SAMPLENAME.$gene.merged.depth)
	frac_covered=$(echo "scale=7; $pos_covered / $total_pos" | bc)

	# dump all data into a file

	printf "$SAMPLENAME\t$gene\t$reads_mapped\t$reads\t$frac_mapped\t$mean_depth\t$pos_covered\t$total_pos\t$frac_covered\n" >> \
	$OUTPUT/$SAMPLENAME\_typing_summary.txt


	#remove all intermediate files
	rm $DIR/typing_mapped/$SAMPLENAME.$gene\_L00*.mapped.sam
	rm $DIR/typing_mapped/$SAMPLENAME.$gene\_L00*.mapped.bam
	rm $DIR/typing_mapped/$SAMPLENAME.$gene\_L00*.mapped.sorted.bam
	rm $DIR/typing_mapped/$SAMPLENAME.$gene.lanebams.txt
	rm $DIR/typing_merged/$SAMPLENAME.$gene.merged.bam 
	rm $DIR/typing_merged/$SAMPLENAME.$gene.merged.depth

done


echo "typing finished"
date
