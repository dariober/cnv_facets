[![Build Status](https://travis-ci.com/dariober/cnv_facets.svg?branch=master)](https://travis-ci.com/dariober/cnv_facets)
[![Language](https://img.shields.io/badge/language-R-brightgreen.svg)](https://cran.r-project.org/)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/dariober/cnv_facets)

Somatic copy variant caller for next generation sequencing data based on the
FACETS package


<!-- vim-markdown-toc GFM -->

* [Installation](#installation)
* [Usage](#usage)
* [Input](#input)
* [Output](#output)
* [References](#references)

<!-- vim-markdown-toc -->

**cnv_facets.R** detects somatic Copy Number Variants (CNVs) from next
generation sequencing data such as whole genome, whole exome and targeted
sequencing experiments. In addition, it estimates tumour purity and ploidy.
cnv_facets.R wraps the [FACETS](https://github.com/mskcc/facets) R package in a
convenient script that can be executed from the command line. 

Installation
============

```
bash setup.sh
```

Usage
=====

```
cnv_facets.R -t tumour.bam -n normal.bam -p cnv_calls.
```

Input
=====

Output
======

References
==========
