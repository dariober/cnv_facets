#
# Copyright (C) 2018 University of Glasgow
#
# Author: Dario Beraldi <dario.beraldi@glasgow.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

library(testthat)

context("Test functions")

source('../bin/cnv_facets.R')

mkdir<- function(){
    unlink('tmp_testthat', recursive= TRUE)
    dir.create('tmp_testthat')
}

teardown({
    unlink('tmp_testthat', recursive= TRUE)
})

test_that("Can write header", {
    header<- make_header(gbuild= 'hg38', genomes= genomes, is_chrom_prefixed= TRUE, cmd= '##Command foo bar', extra= c(purity= 0.5, ploidy= 3.1))
    header<- unlist(strsplit(header, split= '\n'))
    expect_equal(49, length(header))
    expect_true(all(grepl('^##|^#CHROM\t', header)))
    expect_true('##purity=0.5' %in% header)
    expect_true('##ploidy=3.1' %in% header)
    expect_true('##Command foo bar' %in% header)
    contigs<- grep('^##contig=<ID=chr', header, value= TRUE)
    expect_equal(25, length(contigs))
    expect_error(make_header(gbuild= 'xx38', genomes= genomes, is_chrom_prefixed= TRUE, cmd= '##Command foo bar', extra= c(purity= 0.5, ploidy= 3.1)))
})

test_that("Can estimate insert size", {
    avg<- avg_insert_size('data/TCRBOA6-N-WEX.sample.bam')
    expect_equal(168.0, round(avg, 1))
    
    avg<- avg_insert_size('data/single_end.bam', default= 999)
    expect_equal(999, avg)
})

test_that("Can execute pileup on chrom", {
    mkdir()
    exec_snp_pileup(21, 'data/common.sample.vcf.gz', 'tmp_testthat/tmp.cvs.gz', 
        normal_bam= 'data/TCRBOA6-N-WEX.sample.bam', 
        tumour_bam= 'data/TCRBOA6-T-WEX.sample.bam', 
        mapq= 10, 
        baq= 10, 
        pseudo_snp= 250)
    expect_true(file.exists('tmp_testthat/tmp.cvs.gz'))
    expect_equal(0, system2(c('gzip', '--test', 'tmp_testthat/tmp.cvs.gz')))

    csv<- read.table('tmp_testthat/tmp.cvs.gz', sep= ',', header= TRUE)
    expect_equal(12, ncol(csv))
    expect_true(nrow(csv) > 1000)
    expect_equal(21, unique(csv$Chromosome)) 
})

test_that("Pileup thorws error if bash script throws error", {
    mkdir()
    expect_error(
        exec_snp_pileup(21, 'data/common.sample.vcf.gz', 'tmp_testthat/tmp.cvs.gz', 
            normal_bam= 'data/TCRBOA6-N-WEX.sample.bam.bai', 
            tumour_bam= 'data/TCRBOA6-T-WEX.sample.bam', 
            mapq= 10, 
            baq= 10, 
            pseudo_snp= 250)
    )
})

test_that("Can concatenate CSVs", {
    mkdir()
    concat_csv(c('data/stomach.csv.gz', 'data/stomach.csv.gz'), 'tmp_testthat/out.csv.gz', 'tmp_testthat')

    expect_true(file.exists('tmp_testthat/out.csv.gz'))
    expect_equal(0, system2(c('gzip', '--test', 'tmp_testthat/out.csv.gz')))
    csv<- read.table('tmp_testthat/out.csv.gz', sep= ',', header= TRUE)
    expect_equal(599694, nrow(csv))
    expect_equal('Chromosome', names(csv)[1])
    expect_false('Chromosome' %in% csv$Chromosome)
})
