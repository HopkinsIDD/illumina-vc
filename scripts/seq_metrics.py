#!/usr/bin/env python

"""Functions to check read and mapping qulity for a sample based on:
    1. Raw fastq file
    2. unmerged bam file
    3. merged bam file with duplicates marked

"""

import numpy as np
import pandas as pd
import sys
import os
import pysam

from Bio import SeqIO


NAME = sys.argv[1] # sample name

DIR = sys.argv[2] # This must be the same project directory as in other scripts
DIR = DIR.rstrip('/') # trim '/' if necessary 

FASTQ_DIR = sys.argv[3] # absolute path of directory containing raw sequencing data for all samples
FASTQ_DIR = FASTQ_DIR.rstrip('/') # trim '/' if necessary 

NUM_LANE = int(sys.argv[4]) # number of lanes

OUTPUT_DIR = sys.argv[5] # absolute path of directory for reading metrics flie
OUTPUT_DIR = OUTPUT_DIR.rstrip('/') # trim '/' if necessary 

# make directory for output folder if not existing
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)



# initialize a dataframe in which to store metrics data
##
# name: Sample name
# reads_filesize: the total size of all fastq files for this sample
# reads_num: total number of fastq reads
# reads_mapped: the number of fastq reads mapped to the reference genome
# reads_unmapped: the number of fastq reads that were not mapped to the reference genome
# reads_mapped_rmdup: the number of non-duplicated fastq reads mapped to the reference genome
# ref_mean_coverage: mean of the number of reads that covered each position across the genomes
# ref_median_coverage: median of the number of reads that covered each position across the genomes
# ref_pct_coverage: percentage of the genomes that were covered by reads (percentage of non-'N' positions in final assembly)

df = pd.DataFrame(columns=["name","reads_filesize", \
"reads_num","reads_mapped","reads_mapped_rmdup","reads_unmapped", \
"ref_mean_coverage","ref_median_coverage","ref_pct_coverage"])



## file paths
fastq_path = [FASTQ_DIR + "/" + NAME + "_L00" + str(i) + "_R" + str(j) + "_001.fastq.gz" for i in range(1,NUM_LANE+1) for j in range(1,3)]

##CAUTION: The bam file naming convention here followed those in 'align_reference.sh', those should be customized if necessary
bam_unmerged_path = [DIR + "/01_mapped/" + NAME + "_L00" + str(i) + ".mapped.bam" for i in range(1,NUM_LANE+1)]
bam_merged_path = DIR + "/02_merged/" + NAME + ".merged.mkdup.bam"
depth_path = DIR + "/02_merged/" + NAME + ".depth"
genome_path = DIR + "/04_assembly/" + NAME + ".bcftools.fasta"


# get filesize of unfiltered fastq file
reads_filesize = sum([os.path.getsize(fq) for fq in fastq_path])

# get mapping metrics from pysam
reads_unmapped =  sum([int(pysam.view("-c","-f","4",bam_unmerged).rstrip()) for bam_unmerged in bam_unmerged_path]) 
reads_mapped = int(pysam.view("-c","-F","4",bam_merged_path).rstrip())
reads_mapped_rmdup = int(pysam.view("-c","-F","1028",bam_merged_path).rstrip())
reads_num = reads_mapped + reads_unmapped


# get mean and median coverage depth
ref_depth = pd.read_csv(depth_path,header=None,sep='\t',names=["ctg","pos","depth"])
ref_mean_coverage = ref_depth["depth"].mean()
ref_median_coverage = ref_depth["depth"].median()

# get coverage percentage (percent of non-'N' in assembled genomes)
ambigbases = 0
genomelen = 0

with open(genome_path) as genome:
    for record in SeqIO.parse(genome, "fasta"):
        ambigbases+=record.seq.upper().count('N')
        genomelen+=len(record.seq)

ref_pct_coverage = (1 - float(ambigbases*1.0/genomelen)) * 100



# append to dataframe
df.loc[len(df)] = [NAME,reads_filesize,reads_num,reads_mapped,reads_mapped_rmdup,reads_unmapped,\
ref_mean_coverage,ref_median_coverage,ref_pct_coverage]

outfile = OUTPUT_DIR + "/" + NAME + "_testmetrics.csv"
df.to_csv(outfile,index = False)
