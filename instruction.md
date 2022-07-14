## Vibrio cholerae Genome Assembly from Illumina data

This is a guide for taking reads generated on an Illumina sequencer and assembly by aligning reads to the N16961 reference.

## Table of Contents

* [Setting up your working directory](#setting-up-your-working-directory)
* [Installation](#installation)
* [Alignment](#alignment)
* [Variants calling and assembly](#variants-calling-and-assembly)
* [Calculation of run metrics](#calculation-of-run-metrics)
* [Masking](#masking)
* [Illumina contigs alignment](#illumina-contigs-alignment)

## Setting up your working directory

Before you start, you will need to clone this repository. Cloning the repository will ensure you have all the necessary scripts and folders to run the analysis.

First, let's create a project directory to store all the analyses for this sequencing run. This will be your **working directory** for the rest of this process. Open up a terminal window and type the following command:

```
mkdir my-project
```
Where **my-project** can be any name you chose to represent this sequencing run.

To enter this directory, type:

```
cd my-project
```

To get the full path to this directory, now type:

```
pwd
```

The output, which might look something like `/home/username/my-project/`, will be your **working directory** for this analysis. Anytime you are asked to provide the path to your working directory, copy in this information.

From inside your **working directory**, type the following commands to copy this github repository onto your computer:

```
apt update
apt install git
git clone https://github.com/HopkinsIDD/illumina-vc.git .
```

### Finding your raw data

Some of the commands below require you to know where your raw Illumina data is stored. When you started your Illumina run, you selected a data folder in which to store the results. Figure out the path to your raw data and provide it below whenever the command calls for the **path-to-raw-data**.


### Determining the number of processors on your computer

It can be helpful to know how many processors are available on your computer, since using more processors can make some scripts run faster. Run this command to figure out how many you have available:

```
cat /proc/cpuinfo | grep processor
```
Use the bottom number on your screen as the number you provide any time a script askes you for **num-cores**.


## Installation

Before starting, make sure you are in your **working directory**:

```
cd my-project
```

If you have not installed the software needed for assembly, run the following installation script:

```
sudo bash scripts/install-software.sh path-to-software
```

If some of the software is already installed on your computer, make sure that the `path-to-software` matches where the software is currently installed. If this is the first time you are installing assembly software, we reccomend the following path:

```
~/software
```
Please note that you will need an internet connection to complete this step, and that downloading the software will take some time (approximately 30-60min, depending on the speed of your internet connection). If the script is stopped for any reason, you can safely restart it with the command above and it will continue where it left off. If you are ever asked to respond to a yes/no questions, type "y" to continue.


## Alignment

Before starting, make sure you are in your **working directory**:

```
cd my-project
```

First, we need to align the raw Illumina reads to the reference genome running the 'align_reference.sh' script as below. The output of this step is a sorted bam file (bam files from multiple-lane run will be merged into one).

```
bash scripts/align_reference.sh sample-name my-project path-to-raw-data num-lanes num-cores ref-genome picard-path bam-flagstats
```

To determine the path to enter instead of path-to-raw-data, see [Finding your raw data](#finding-your-raw-data) above. **This script assumes the inputs are pair-ended Illumina reads with standard [Illumina naming convention](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm). _If your fastq files were not in this convention, please rename them first._** The `sample-name` is the `SampleName` indicated in the link. `num-lanes` is the number of lanes you used to run this sample. `num-cores` should be less or equal to [the number of processors on your computer](#determining-the-number-of-processors-on-your-computer). The path to `ref-genome` is "ref_genomes/vc_reference.fasta". If you installed the picard software using the 'install-software.sh' script, the `picard-path` is "path-to-software/picard/build/libs/picard.jar". If you have it installed yet in your computer, check the path by:

```
find path-to-software -name picard.jar
```

Where `path-to-software` is the path of your installed softwares.

Lastly, if you want to calculate the alignment statistics for the bam file, you can set `bam-flagstats` to be 1, otherwise 0.

## Variants calling and assembly

Before starting, make sure you are in your **working directory**:

```
cd my-project
```

We are going to identify the variants from the reference genome and create consensus genome for this sample. Run the script as below:

```
bash scripts/call_variants.sh sample-name my-project path-to-input-bam num-cores ref-genome snp-support req-depth req-strand-depth
```

For `sample-name`, `my-project`, `num-scores`, `ref-genome`, check instructions in [Alignment](#alignment). `path-to-input-bam` is the path to the input bam file. If you run this script following 'align_reference.sh', `path-to-input-bam` is '`my-project`/02_merged/`sample-name`.merged.mkdup.bam'. 

**`snp-support, req-depth and req-strand-depth`** are variant calling parameters. Their definitions were listed below:

**`snp-support`**: The percent(0-1) of mapped reads required to support alternate allele.

**`req-depth`**: The minimum total depth required for each base.

**`req-strand-depth`**: The minumum number of reads matching alternate allele per strand.

## Calculation of run metrics

Before starting, make sure you are in your **working directory**:

```
cd my-project
```

After running [Alignment](#alignment) and [Variant calling](#variants-calling-and-assembly) scripts, you can calculate some run metrics as below:

```
python scripts/seq_metrics.py sample-name my-project path-to-raw-data num-lanes out-dir
```
For `sample-name`, `my-project`, `path-to-raw-data`, `num-lanes`, check instructions in [Alignment](#alignment). `out-dir` is user-defined metrics results output directory. Below shows the metrics calculated by this script:

**name**: Sample name

**reads_filesize**: the total size of all fastq files for this sample

**reads_num**: total number of fastq reads

**reads_mapped**: the number of fastq reads mapped to the reference genome

**reads_unmapped**: the number of fastq reads that were not mapped to the reference genome

**reads_mapped_rmdup**: the number of non-duplicated fastq reads mapped to the reference genome

**ref_mean_coverage**: mean of the number of reads that covered each position across the genomes

**ref_median_coverage**: median of the number of reads that covered each position across the genomes

**ref_pct_coverage**: percentage of the genomes that were covered by reads (percentage of non-'N' positions in final assembly)


## Masking

Before starting, make sure you are in your **working directory**:

```
cd my-project
```

In this step, we will use a previously-published list of potential recombinant sites to mask areas of the genome that could mess up our phylogenetic analysis. The `recomb_mask.gff ` file used to mask sites was originally published at https://figshare.com/s/d6c1c6f02eac0c9c871e by Weill et al. We can run the masking script with the following command:

```
bash scripts/mask_fasta.sh sample-name my-project path-to-consensus-genome recomb_mask.gff  chrom-name
```

For `sample-name`, `my-project`, check instructions in [Alignment](#alignment). `path-to-consensus-genome` is the path of the consensus genome fasta files. If you run this script following 'call_variants.sh', `path-to-consensus-genome` is '`my-project`/04_assembly/`sample-name`.bcftools.fasta'. `chrom-name` is the chromosome name in the gff file. It's 'gff_seqname' in this file. **You can check that in the first column of other gff file**.


These genomes now have recombinant regions masked, and can be used in phylogenetic analyses. For best results, run gubbins on the final, masked alignment prior to phylogenetic analysis. 

## Illumina contigs alignment

Sometimes, only contigs rather than raw Illumina reads are available to us. We are going to run a software 'snippy' to align contigs to the reference genome and obtain the bam file for downstream processing:

```
bash scripts/run_snippy sample-name my-project path-to-contigs snippy-output num-cores ref-genome
```

`sample-name` is defined by user. For  `my-project`, `num-scores`, `ref-genome`, check instructions in [Alignment](#alignment). `path-to-contigs` is the path to contigs file. `snippy-output` is the snippy output folder.

The path to the output bam file is '`snippy-output`/`sample-name`/snps.bam'. Then you can put this in [Variant calling](#variants-calling-and-assembly). **Warning: [Run metrics](#calculation-of-run-metrics) can not be run with contigs file.** Also, we suggest you create a separate working directory for contigs assembly to avoid messing files up.