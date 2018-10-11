#!/usr/bin/env bash

set -e 
set -o pipefail

rm -rf test_out
../bin/cnv_facets.R -p data/pileup.csv.gz -o test_out/out

ls test_out/out.vcf.gz
ls test_out/out.cnv.png
ls test_out/out.spider.pdf

rm -rf test_out
../bin/cnv_facets.R -p data/pileup.csv -o test_out/out

ls test_out/out.vcf.gz
ls test_out/out.cnv.png
ls test_out/out.spider.pdf

rm -rf test_out
../bin/cnv_facets.R -t data/tumour.bam -n data/normal.bam -vcf data/snps.vcf.gz -o test_out/out

ls test_out/out.vcf.gz
ls test_out/out.cnv.png
ls test_out/out.spider.pdf
ls test_out/out.csv.gz

rm -rf test_out
../bin/cnv_facets.R -p data/stomach.csv.gz -o test_out/out

ls test_out/out.vcf.gz
ls test_out/out.cnv.png
ls test_out/out.spider.pdf
ls test_out/out.csv.gz

rm -rf test_out

