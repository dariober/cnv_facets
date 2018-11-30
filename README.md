[![Build Status](https://travis-ci.org/wwcrc/cnv_facets.svg?branch=master)](https://travis-ci.org/wwcrc/cnv_facets)
[![Language](https://img.shields.io/badge/language-R-brightgreen.svg)](https://cran.r-project.org/)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/wwcrc/cnv_facets)

Somatic copy variant caller for next generation sequencing data based on the
[FACETS](https://github.com/mskcc/facets) package

<!-- vim-markdown-toc GFM -->

* [Requirements and Installation](#requirements-and-installation)
* [Input](#input)
    * [Option 1: BAM & VCF input](#option-1-bam--vcf-input)
    * [Option 2: Pileup input](#option-2-pileup-input)
* [Output](#output)
    * [Variants](#variants)
    * [CNV profile plot](#cnv-profile-plot)
    * [Diagnostic plot](#diagnostic-plot)
    * [Pileup file](#pileup-file)
* [Usage guidelines](#usage-guidelines)
* [Time & memory footprint](#time--memory-footprint)

<!-- vim-markdown-toc -->

**cnv_facets.R** detects somatic Copy Number Variants (CNVs) from next
generation sequencing data such as *whole genome*, *whole exome* and *targeted
sequencing* experiments. In addition, it estimates tumour purity and ploidy.
cnv_facets.R is based on the original package [FACETS](https://github.com/mskcc/facets)

The advantage of **cnv_facets.R** over the original [FACETS](https://github.com/mskcc/facets) 
package is in the convenience of executing all the necessary steps, from BAM input to VCF 
output, in a single command line call.

Requirements and Installation
=============================

cnv_facets.R requires a reasonably recent version of
[R](https://cran.r-project.org/) on a `*Nix` operating system. At the time of
this writing, it has been developed and deployed on R 3.5 on CentOS 7.

To compile and install execute:

```
bash setup.sh --bin_dir </dir/on/path>
```

Where `/dir/on/path` is a directory on your PATH where you have permission to
write, *e.g.*, `~/bin`.

`setup.sh` accomplishes three main tasks:

* Install any missing, required R package

* Compile the helper program `snp-pileup` and move it to the designated PATH
  directory unless already found on PATH

* Run the test suite

Input
=====

Option 1: BAM & VCF input
-------------------------

* A bam file of the **tumour** sample

* A bam file of the **normal** sample (matched to the tumour, typically a blood
  sample from the same patient)

* A VCF file of common, polymorphic SNPs. For human samples, a good source is
  the dbSNP file
  [common_all.vcf.gz](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/). 
  See also NCBI [human variation sets in VCF Format](https://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/).

**USAGE**

```
cnv_facets.R -t <tumour.bam> -n <normal.bam> -vcf <snps.vcf.gz> -o <output_prefix> [...]
```

Option 2: Pileup input
----------------------

* A pileup file, comma separated, of read counts for the reference and alternate allele
  at polymorphic SNPs. This file must have the following columns (order of
  columns is not important, additional columns are ignored):
 
 * *Chromosome* Chromosome of the SNP
 
 * *Position* Position of the SNP
 
 * *File1R* Read depth supporting the REF allele in **normal** sample
 
 * *File1A* Read depth supporting the ALT allele in **normal** sample
 
 * *File2R* Read depth supporting the REF allele in **tumour** sample
 
 * *File2A* Read depth supporting the ALT allele in **tumour** sample

These are the first lines of the test file `test/data/stomach.csv.gz`
accompanying the original facets package:

```
"Chromosome","Position","Ref","Alt","File1R","File1A","File1E","File1D","File2R","File2A","File2E","File2D"
1,69424,N,N,170,117,0,0,158,103,0,0
1,69515,N,N,0,76,0,0,0,77,0,0
1,69536,N,N,103,0,0,0,99,0,0,0
1,808866,N,N,96,0,0,0,133,0,0,0
1,809120,N,N,66,0,0,0,105,0,0,0
```

This pileup file is generated by `cnv_facets.R` when run with bam input as in
option 1. Alternatively, it can be created by running `snp-pileup`. Using the
pileup file as input instead of the bam files has the advantage of saving
computing time.

**USAGE**

```
cnv_facets.R -p <pileup.csv.gz> -o <output_prefix> [...]
```

Output
======

The option `--out/-o <prefix>` determines the name and location of the output
files. For more information refer to the documentation of the FACETS package.

Variants
--------

* `<prefix>.vcf.gz`

VCF file compressed and indexed of copy number variants. The INFO tags below annotate each variant:

Tag | Type | Description
----|------|------------
SVTYPE | String | Type of structural variant
SVLEN | Integer | Difference in length between REF and ALT alleles
END | Integer | End position of the variant described in this record
NUM_MARK | Integer | Number of SNPs in the segment
NHET | Integer | Number of SNPs that are deemed heterozygous
CNLR_MEDIAN | Float | Median log-ratio (logR) of the segment. logR is defined by the log-ratio of total read depth in the tumor versus that in the normal
CNLR_MEDIAN_CLUST | Float | Median log-ratio (logR) of the segment cluster. logR is defined by the log-ratio of total read depth in the tumor versus that in the normal
MAF_R | Float | Log-odds-ratio (logOR) summary for the segment. logOR is defined by the log-odds ratio of the variant allele count in the tumor versus in the normal
MAF_R_CLUST | Float | Log-odds-ratio (logOR) summary for the segment cluster. logOR is defined by the log-odds ratio of the variant allele count in the tumor versus that in the normal
SEGCLUST | Integer | Segment cluster to which the segment belongs
CF_EM | Float | Cellular fraction, fraction of DNA associated with the aberrant genotype. Set to 1 for normal diploid
TCN_EM | Integer | Total copy number. 2 for normal diploid
LCN_EM | Integer | Lesser (minor) copy number. 1 for normal diploid
CNV_ANN | String | Annotation features assigned to this CNV

The header of the VCF file also stores the estimates of tumour purity and ploidy.

CNV profile plot
----------------

* `<prefix>.cnv.png`

Summary plot of CNVs across the genome, for [example](./docs/tex.cnv.png):

<img src="./docs/tex.cnv.png" height="600"/>

Diagnostic plot
---------------

* `<prefix>.spider.pdf`

This is a diagnostic plot to check how well the copy number fits
work The estimated segment summaries are plotted as circles
where the size of the circle increases with the number of loci in
the segment. The expected value for various integer copy number
states are drawn as curves for purity ranging from 0 to 0.95. For
a good fit, the segment summaries should be close to one of the
lines. (*Description from `facets::logRlogORspider`*). For [example](./docs/tex.spider.png):

<img src="./docs/tex.spider.png" height="400"/>

Pileup file
-----------

* `<prefix>.csv.gz`

File of nucleotide counts at each SNP in normal and tumour sample.

Usage guidelines
================

* Option `--cval`

Critical values for segmentation in pre-processing and processing.
Larger values reduce segmentation. [25 150] is facets default based on exome data. For whole genome
consider increasing to [25 400] and for targeted sequencing consider reducing them. Default 25 150

* Option `--nbhd-snp`

If an interval of size nbhd-snp contains more than one SNP, sample a random one.
This sampling reduces the SNP serial correlation. This value should be similar
to the median insert size of the libraries. 250 is facets default based on
exome data. For whole genome consider increasing to 500 and for target
sequencing decrease to 150. Default 250

Time & memory footprint
=======================

The analysis of a whole genome sequence where the
tumour is sequenced at ~80x (~2 billion reads, BAM file ~200 GB) and the normal
at ~40x (~1 billion reads, BAM files ~100 GB) with ~37 million SNPs (from dbSNP
`common_all_20180418.vcf.gz`) and with no filtering on read depth and read
quality requires:

* 5 hours to prepare the SNP pileup with small memory footprint. Time is mostly
  driven by the size of the BAM files. To speed-up the pileup consider the 
  option `--ncores` to parallelize across chromosomes.

* 1 hour and ~15 GB of memory for the actual detection of CNVs starting from
  the pileup. Time and memory is mostly driven by the number of SNPs

This is A typical targeted sequencing datasets wi

