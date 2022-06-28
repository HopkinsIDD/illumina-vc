#!/bin/bash

echo "Beginning of script"
date

#use snippy to obtain alignment from contigs
#In the output directory, the file named 'snps.bam' could be input into 'call_variant.sh'


# ---------------- Input Paramters and Files  ---------------- #

# path to project directory 
DIR=$1
DIR=${DIR%/} # remove trailing slash from directory name if necessary

# path to contigs file
CONTIG=$2

# path to snippy output directory
OUTPUT=$3
OUTPUT=${OUTPUT%/}

# sample name 
SAMPLENAME=$4

# number of available cores
NUM_CORES=$5

# path to reference genome
REF_GENOME=$6


cd $DIR
echo $SAMPLENAME

# create folder to store tmp files for checking file completion

if [ ! -d check_tmp_snippy ]; then
	mkdir check_tmp_snippy
fi


printf "running snippy for: %s" "$SAMPLENAME"



# ---------------- Running snippy  ---------------- #

date

#check if snippy has been run on this sample
if [ ! -d $OUTPUT/$SAMPLENAME ]; then
	snippy --outdir $OUTPUT/$SAMPLENAME --ctgs $CONTIG --ref $REF_GENOME --cpus $NUM_CORES
	touch check_tmp_snippy/.$SAMPLENAME.snippy.check
elif [ ! -e check_tmp_snippy/.$SAMPLENAME.snippy.check ]; then
	rm -r $OUTPUT/$SAMPLENAME
	snippy --outdir $OUTPUT/$SAMPLENAME --ctgs $CONTIG --ref $REF_GENOME --cpus $NUM_CORES
	touch check_tmp_snippy/.$SAMPLENAME.snippy.check	
fi

echo "snippy run completed"

date
