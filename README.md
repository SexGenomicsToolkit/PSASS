[![GitHub release (latest by date)](https://img.shields.io/github/v/release/RomainFeron/PSASS?color=lightorange)](https://github.com/RomainFeron/PSASS/releases)
[![Conda (channel only)](https://img.shields.io/conda/vn/bioconda/psass?color=lightorange)](https://bioconda.github.io/recipes/psass/README.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3702337.svg)](https://doi.org/10.5281/zenodo.3702337)

# Pooled Sequencing Analyses for Sex Signal (PSASS)

Current release: 3.2.0

## Overview

PSASS (Pooled Sequencing Analysis for Sex Signal) is a software to compare pooled sequencing datasets from two groups (usually two sexes). Results from PSASS can be easily visualized using the [sgtr](https://github.com/SexGenomicsToolkit/sgtr) R package. PSASS is integrated in a [Snakemake workflow](https://github.com/SexGenomicsToolkit/PSASS-workflow) to perform all required steps starting from a genome and reads files. PSASS was developed as part of a project by the [LPGP](https://www6.rennes.inra.fr/lpgp/) lab from INRA, Rennes, France.

**Citing PSASS**

There is currently no paper officially describing PSASS. To cite PSASS, use the DOI provided by [Zenodo](https://doi.org/10.5281/zenodo.3702337).

## Installation

### Install with conda

PSASS is available in [Bioconda](https://bioconda.github.io/recipes/psass/README.html). To install psass with conda, run the following command:

```bash
conda install -c bioconda psass
```

### Install from source

PSASS implements parsing of alignment files using the [htslib](https://github.com/samtools/htslib) library, which requires **autotools** to build and depends on **zlib**, **libbz2**, and **liblzma** to read CRAM files. All are installed by default on most linux distributions. Compilation was tested with **gcc >= 5.3.0**.

To build psass, follow these instructions:

```bash
# Clone the repository
git clone https://github.com/RomainFeron/PSASS.git
# Navigate to the PSASS directory
cd psass
# Build PSASS
make
```

## Usage

In summary, *psass* takes as input one reads alignment file for each pool and computes the following metrics:

- Between-pools Fst in a sliding window
- Number of SNPs specific to each pool, defined as SNPs heterozygous in one pool and homozygous in the other, in a sliding window
- Absolute and relative depth for each pool in a sliding window
- (optional) The position of all bases with high between-pools Fst
- (optional) The position of all SNPs specific to each pool
- (optional) Number of SNPs specific to each pool, and absolute and relative depth in each pool for all genes in a provided GFF file, for both coding and noncoding parts of the genes

Currently, *psass* implements three commands:

- `pileup` : generate a nucleotides depth file from two alignment files
- `analyze` : compute pool comparison metrics from a nucleotides depth file
- `convert` : convert output from samtools mpileup to a nucleotides depth file (deprecated, use `pileup` instead)


### Quickstart

In this example, we will compare a pool of female individuals and a pool of male individuals with `psass`. We assume the following input data:

- *genome.fa*: the assembly to which reads from each pool were aligned
- *females.cram*: alignment file for the reads from the female pool (sorted by genomic coordinates)
- *males.cram*: alignment file for the reads from the male pool (sorted by genomic coordinates)

#### Generate a pileup file with `psass pileup`:

```bash
psass pileup --reference genome.fa --output-file pileup.tsv females.cram males.cram
```

This command generates the file *pileup.tsv* which contains the nucleotide composition of each base in *genome.fa* for each pool.

#### Compute metrics with `psass analyze`:

```bash
psass analyze --window-size 10000 --output-resolution 1000 --snp-file psass_snps.tsv pileup.tsv psass_window.tsv
```

This command generates two output files:
- *psass_window.tsv* which contains between-pool FST, pool-specific SNPs, and depth for each pool in a sliding window of 10,000 bp, output every 1,000 bp.
- *psass_snps.tsv* which contains the position of each pool-specific SNP.


### pileup

The **pileup** function generates a file with the nucleotide composition of all genomic positions for any number of alignment files. Alignment files can be provided in CRAM or BAM format and need to be sorted by coordinate. The output is a wig-like file with a header line giving the contig name and length followed by one line giving the nucleotide composition in each alignment file for each positions in the contig for all contigs in the alignment files.

```
Usage: psass pileup [OPTIONS] ALIGNMENT_FILES...
```

**Arguments**:

Argument               | Type       | Description                                                      | Default |
-----------------------|------------|------------------------------------------------------------------|---------|
ALIGNMENT_FILES        |  `file`    |  One alignment file for each pool, in either CRAM or BAM format  |         |
--reference, -r        |  `file`    |  Reference file in fasta format, required with CRAM input files  |         |
--output-file, -o      |  `string`  |  Write output to this file instead of stdout                     |         |
--min-map-quality, -q  |  `string`  |  Minimum mapping quality to include a read in pileup             |    0    |
--help                 |            |  Display help message                                            |         |

**pileup** generates a file with the following format:

```
#Files   <first_input_file>   <second_input_file>   # Comment line
region=<contig name>   len=<contig length>          # Header line for the first contig, encoding the contig name and its length
nA,nT,nC,nG,nN,nO   nA,nT,nC,nG,nN,nO               # Count for each type of nucleotide (comma-separated) in each pool (tab-separated) for pos 0
nA,nT,nC,nG,nN,nO   nA,nT,nC,nG,nN,nO               # Count for each type of nucleotide (comma-separated) in each pool (tab-separated) for pos 1
...
region=<contig name>   len=<contig length>          # Header line for the second contig, encoding the contig name and its length
nA,nT,nC,nG,nN,nO   nA,nT,nC,nG,nN,nO               # Count for each type of nucleotide (comma-separated) in each pool (tab-separated) for pos 0
nA,nT,nC,nG,nN,nO   nA,nT,nC,nG,nN,nO               # Count for each type of nucleotide (comma-separated) in each pool (tab-separated) for pos 1
...
```

### analyze

The **analyze** function computes between-pool FST, pool-specific SNPs, and depth in a sliding window from a nucleotide composition file generate with `psass pileup` from two alignment files.

```
Usage: psass analyze [OPTIONS] INPUT_FILE OUTPUT_FILE
```

**Arguments**:

Argument                 | Type       | Description                                                             |  Default  |
-------------------------|------------|-------------------------------------------------------------------------|-----------|
INPUT_FILE               |  `file`    |  Path to a nucleotides depth file generated by psass pileup or convert  |           |
OUTPUT_FILE              |  `string`  |  Path to an output file for sliding window metrics                      |           |
--pool1, -p              |  `string`  |  Name of the first pool                                                 |  females  |
--pool2, -q              |  `string`  |  Name of the second pool                                                |  males    |
--snp-file, -s           |  `string`  |  Output pool-specific SNPs to this file                                 |           |
--fst-file, -f           |  `string`  |  Output high FST positions to this file                                 |           |
--genes-file, -g         |  `string`  |  Output gene metrics to this file (requires a GFF file)                 |           |
--gff-file, -G           |  `string`  |  Path to a GFF file for gene-specific output                            |           |
--popoolation            |            |  If set, assumes the input file was generated with popoolation2         |           |
--min-depth, -d          |  `int`     |  Minimum depth to include a site in the analyses                        |   10      |
--window-size, -w        |  `int`     |  Size of the sliding window (in bp)                                     |   100000  |
--output-resolution, -r  |  `int`     |  Output resolution for sliding window metrics (in bp)                   |   10000   |
--freq-het, -e           |  `float`   |  Allele frequency to consider a SNP heterozygous in a pool              |   0.5     |
--range-het, -u          |  `float`   |  Range of allele frequency to consider a SNP heterozygous in a pool     |   0.15    |
--freq-hom, -o           |  `float`   |  Allele frequency to consider a SNP homozygous in a pool                |   1       |
--range-hom, -v          |  `float`   |  Range of allele frequency to consider a SNP homozygous in a pool       |   0.05    |
--min-fst, -t            |  `float`   |  Minimum FST to output a site in the FST positions file                 |   0.1     |
--group-snps             |            |  If set, group consecutive snps to count them as a single polymorphism  |           |
--help                   |            |  Print this help message and exit                                       |           |

#### Output files

**Sliding window output file**

A tabulated file with contig, position on the contig, contig length, number of pool-specific SNPs, between-pool FST in the window, absolute depth, and relative depth for each pool in a sliding window of size given by `--window-size`. Output every *N* bp, with *N* given by `--output-resolution`.

```
Contig   Position  Length  Snps_<pool1>  Snps_<pool2>       Fst  Abs_depth_<pool1>  Abs_depth_<pool2>  Rel_depth_<pool1>  Rel_depth_<pool2>
Contig1         0    6000             4             5    0.0000                166                174               0.73               0.74
Contig1     10000    6000             4             6    0.0000                156                165               0.68               0.70
Contig1     20000    6000             5             6    0.0000                181                193               0.79               0.82
Contig1     30000    6000             5             6    0.0000                167                178               0.74               0.76
Contig1     40000    6000             6             6    0.0000                154                164               0.68               0.70
...
```

**Pool-specific SNPs output file**

A tabulated file with contig, position on the contig, pool ID, and nucleotides frequencies at this position in each pool (O for 'other') for all positions called as pool-specific SNPs.

```
Contig   Position  Pool     <pool1>_A  <pool1>_T  <pool1>_C  <pool1>_G  <pool1>_N  <pool1>_O  <pool2>_A  <pool2>_T  <pool2>_C  <pool2>_G  <pool2>_N  <pool2>_I
Contig1      6170  <pool1>       0.06       0.19       0.12       0.62       0.00       0.00       0.05       0.00       0.00       0.95       0.00       0.00
Contig1      6176  <pool2>       0.62       0.08       0.23       0.08       0.00       0.00       1.00       0.00       0.00       0.00       0.00       0.00
Contig1      7090  <pool2>       1.00       0.00       0.00       0.00       0.00       0.00       0.65       0.18       0.12       0.06       0.00       0.00
Contig1      7101  <pool1>       0.00       0.00       0.05       0.95       0.00       0.00       0.04       0.13       0.22       0.61       0.00       0.00
...
```

**High FST positions output file**

A tabulated file with contig, position on the contig, and between-pool FST at this position for all positions with FST higher than the threshold specified with `--min-fst`.

```
Contig   Position     Fst
Contig1      4156  0.3333
Contig1      4157  0.6000
Contig1      4158  0.6000
Contig1      4159  0.3333
...
```

#### Principle

In this section, we will briefly describe the process implemented in PSASS to compute FST and pool-specific SNPs.

The input of PSASS `analyze` is a file that contains nucleotide counts for each pool and each position in the reference genome. PSASS parses this input and computes:

* FST between the two pools for each position, using the formula implemented in [Popoolation2](https://academic.oup.com/bioinformatics/article/27/24/3435/306737).
* FST between pools in a sliding window, using the estimate defined in [Karlsson *et al.*, 2007](https://www.nature.com/articles/ng.2007.10) (in Supplementary information).
* The position of pool-specific SNPs. A pool-specific SNP is called at a given genomic position if, for any possible nucleotide, the frequency of this nucleotide is 1 in one pool and 0.5 in the other pool. Values for these thresholds can be adjusted with command-line arguments, and a range can be given; by default these values are: 1 - [0-0.05] in the homozygous pool, and 0.5 +/- 0.1 in the heterozygous pool).
* Number of pool-specific SNPs in a sliding window, using the definition from above.

Therefore, the definition of FST and pool-specific SNPs used in PSASS is only base on comparing the pools themselves; the reference used for the alignments is not used to call SNPs.

### convert

The **convert** function converts a pileup file generated by samtools mpileup into a psass pileup file. It has been replaced by the `psass pileup` command.

```
Usage: psass convert [OPTIONS] INPUT
```

**Arguments**:

Argument           | Type       | Description                                                         | Default |
-------------------|------------|---------------------------------------------------------------------|---------|
INPUT              |  `file`    |  Either the path to a samtools pileup output file or "-" for stdin  |         |
--output-file, -o  |  `string`  |  Write output to this file instead of stdout                        |         |
--help             |            |  Display help message                                               |         |
