## Vibrio cholerae typing from Illumina data

This is a guide for taking reads generated on an Illumina sequencer and mapping them against specific genes (_ctxA_, _wbeO1_, _wbfO139_, _tcpA_ classical, _tcpA_ El Tor, _toxR_) for the purpose of determining sample serotype. 

## Setting up your working directory

This script assumpes the following:

* That you have already set up your working directory as described in the [assembly instructions](https://github.com/HopkinsIDD/illumina-vc/blob/main/instruction.md). Make sure you know where the raw fastq files are stored and make sure they follow the [naming convention](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm)

* That the `scripts` folder inside your working directory contains the `map_genes.sh` script.

* That all six gene-specific reference FASTAs are contained in a directory called `typing`. If you cloned this git repository, these reference files can be found in `ref_genomes/typing`.

Once you have the required files in your working directory, navigate to this directory and create a new folder where you will store your typing results:

```
cd my-project
```

```
mkdir typing_results
```

## Running the typing script


Run the typing script as follows:

```
bash scripts/map_genes.sh sample-name my-project path-to-raw-data num-lanes num-cores ref-dir output-dir
```

Where `sample-name` is the 'SampleName' prior to '_L00...' in your fastq files,`my-project` is the absolute path to your project directory, and `path-to-raw-data` is the path to your raw fastq files. `num-lanes` is the number of lanes you used to run this sample. `num-cores` should be less or equal to [the number of processors on your computer](https://github.com/HopkinsIDD/illumina-vc/blob/main/instruction.md/#determining-the-number-of-processors-on-your-computer).`ref-dir` should be the directory that contains the gene-specific FASTA files (usually `my-project/ref_genomes/typing`). `output-dir` will be the folder to store typing results (`my-project/typing_results`).




## Evaluating typing results

The typing script will produce one output file (`samplename_typing_summary.txt`). This file will look something like this:

```
SAMPLENAME      gene    reads_mapped    reads   frac_mapped     mean_depth      pos_covered     total_pos       frac_covered
ERR3039943      ctxA    96      588182  .0001632        27.4582 777     777     1.0000000
ERR3039943      tcpA_classical  87      588182  .0001479        9.74046 1048    1048    1.0000000
ERR3039943      tcpA_eltor      96      588182  .0001632        28.3126 675     675     1.0000000
ERR3039943      toxR    134     588182  .0002278        30.8923 891     891     1.0000000
ERR3039943      wbeO1   3660    588185  .0062225        38.4333 21122   21122   1.0000000
ERR3039943      wbfO139 1923    588187  .0032693        9.17084 46721   46721   1.0000000
```

These files may be easier to view in a spreadsheet software such as Excel.
