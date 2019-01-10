#
# Copyright (C) 2018-2019 University of Glasgow
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

test_that("Can annotate intervals", {
    
    # One CNV multiple annotations
    cnv<- data.table(chrom= 'chr1', start= 29649111, end= 29662081)
    annotate(cnv, 'data/annotation.bed')
    expect_equal('C,D', cnv$annotation)
    
    # CNV is neutral
    cnv<- data.table(chrom= 'chr1', start= 29649111, end= 29662081, type= 'NEUTR')
    annotate(cnv, 'data/annotation.bed')
    expect_true(is.na(cnv$annotation))

    # CNV does not overlap any annotation
    cnv<- data.table(chrom= c('chr1', 'chr18'), start= c(1, 1), end= c(100, 100))
    annotate(cnv, 'data/annotation.bed')
    expect_true(all(is.na(cnv$annotation)))

    # One annotation maps to multiple CNVs
    cnv<- data.table(chrom= c('chr2', 'chr2'), start= c(1, 1000), end= c(10, 1010))    
    annotate(cnv, 'data/annotation.bed')
    expect_equal(c('F', 'F'), cnv$annotation)

    # Annotation record has no actual annotation in col 4
    cnv<- data.table(chrom= c('chr4', 'chr4'), start= c(1, 1000), end= c(10, 1010))    
    annotate(cnv, 'data/annotation.bed')
    expect_true(all(is.na(cnv$annotation)))

    cnv<- data.table(chrom= c('chr1', 'chr18', 'chr4'), 
                     start= c(69000, 1, 20000), 
                     end=   c(70000, 10, 20000)
          )
    annotate(cnv, 'data/annotation.bed')
    expect_equal(c('chr1', 'chr18', 'chr4'), cnv$chrom)
    expect_equal(c('A,B,gene%3Db%3BGene%2CFoo', NA, NA), cnv$annotation)
})

test_that("Can run facets", {
    rcmat<- fread('gzip -c -d data/rcmat.txt.gz')
    facets<- run_facets(
           pre_rcmat= rcmat,
           pre_gbuild= 'hg38', 
           pre_snp.nbhd= 250,
           pre_het.thresh= 0.25, 
           pre_cval= 25, 
           pre_deltaCN= 0, 
           pre_unmatched= FALSE, 
           pre_ndepth= 1,
           pre_ndepthmax= 1e8,
           proc_cval= 250, 
           proc_min.nhet= 15, 
           proc_dipLogR= NULL,
           emcncf_unif= FALSE, 
           emcncf_min.nhet= 15,
           emcncf_maxiter= 20,
           emcncf_eps= 1e-3)
    expect_true(is.data.table(facets$emcncf_fit$cncf))
    expect_true(is.null(names(facets$emcncf_fit$purity)))
})

test_that("Can convert facets record to VCF", {
    x<- data.table(seg= 1, 
                   chrom= 'chr1',
                   num.mark=10, 
                   nhet=20, 
                   cnlr.median=0.1, 
                   mafR=0.2, 
                   segclust= 2, 
                   cnlr.median.clust= 0.3, 
                   mafR.clust= 0.4,
                   start= 100,
                   end= 200,
                   cf.em= 0.8,
                   tcn.em= 2,
                   lcn.em= 1,
                   type= 'NEUTR',
                   annotation= NA)
    vcf<- facetsRecordToVcf(x)
    expect_equal(8, length(vcf))
    expect_equal('chr1', vcf[1])
    expect_true(grepl('SVTYPE', vcf[8]))
    expect_true(grepl('CNV_ANN=.', vcf[8]))
    x$annotation<- 'ACTB,KRAS'
    vcf<- facetsRecordToVcf(x)
    expect_true(grepl('CNV_ANN=ACTB,KRAS', vcf[8]))
})

test_that("Can reset chroms", {
    cncf<- data.table(chrom= rep(1:23, each= 2))
    reset_chroms(cncf, gbuild= 'hg38', chr_prefix= FALSE)
    expect_equal(rep(c(1:22, 'X'), each= 2), cncf$chrom) 

    cncf<- data.table(chrom= rep(1:23, each= 2))
    reset_chroms(cncf, gbuild= 'hg38', chr_prefix= TRUE)
    expect_equal(paste0('chr', rep(c(1:22, 'X'), each= 2)), cncf$chrom)

    cncf<- data.table(chrom= rep(1:20, each= 2))
    reset_chroms(cncf, gbuild= 'mm10', chr_prefix= TRUE)
    expect_equal(paste0('chr', rep(c(1:19, 'X'), each= 2)), cncf$chrom)

    cncf<- data.table(chrom= rep(1:21, each= 2))
    expect_error(
        reset_chroms(cncf, gbuild= 'mm10', chr_prefix= TRUE)
    )

    cncf<- data.table(chrom= c(1:22, 'X', 'Y'))
    expect_error(
        reset_chroms(cncf, gbuild= 'hg38', chr_prefix= TRUE)
    )

    cncf<- data.table(chrom= c(1:22, 'X'))
    expect_error(
        reset_chroms(cncf, gbuild= 'FOO', chr_prefix= TRUE)
    )
})

test_that("Can read SNP pileup", {
    rcmat<- readSnpMatrix2('data/stomach.csv.gz', 'hg38')
    expect_false(rcmat$chr_prefix)
    expect_true(all(c('Chromosome', 'Position', 'NOR.DP', 'NOR.RD', 'TUM.DP', 'TUM.RD') == names(rcmat$pileup)))
    expect_equal(299847, nrow(rcmat$pileup))
    
    system('gzip -c -d data/stomach.csv.gz > data/stomach.csv')
    rcmat<- readSnpMatrix2('data/stomach.csv.gz', 'hg38')
    unlink('data/stomach.csv')
    expect_equal(299847, nrow(rcmat$pileup))

    rcmat<- readSnpMatrix2('data/stomach_chr.csv.gz', 'hg38')
    expect_true(rcmat$chr_prefix)
})

test_that("Can filter read count matrix", {
    rcmat<- fread('gzip -c -d data/rcmat.txt.gz')
    flt<- filter_rcmat(rcmat, min_ndepth= 60, max_ndepth= 200, target= NULL)
    expect_equal(60, min(flt$NOR.DP))    
    expect_equal(199, max(flt$NOR.DP))    
    expect_true(nrow(flt) > 1000)
    
    flt<- filter_rcmat(rcmat, min_ndepth= 115, max_ndepth= 2000, target= 'data/targets.bed')
    expect_true(all(c('Chromosome', 'Position', 'NOR.DP', 'NOR.RD', 'TUM.DP', 'TUM.RD') == names(flt)))
    expect_equal(115, min(flt$NOR.DP))
    expect_equal(69515, min(flt$Position))
    expect_equal(879546, max(flt$Position))
    expect_equal(10, nrow(flt))

    ## Test overlap or depth returns 0 positions
    flt<- filter_rcmat(rcmat, min_ndepth= 10000, max_ndepth= 11000, target= 'data/targets.bed')
    expect_true(all(c('Chromosome', 'Position', 'NOR.DP', 'NOR.RD', 'TUM.DP', 'TUM.RD') == names(flt)))
    expect_equal(0, nrow(flt))

    flt<- filter_rcmat(rcmat[Position > 1000000], min_ndepth= 1, max_ndepth= 11000, target= 'data/targets.bed')
    expect_true(all(c('Chromosome', 'Position', 'NOR.DP', 'NOR.RD', 'TUM.DP', 'TUM.RD') == names(flt)))
    expect_equal(0, nrow(flt))
})

test_that("Can classify CNVs", {
    cnv<- fread('data/cnv_classification.tsv')
    classify_cnv(cnv)
    expect_true(all(cnv$cnv_type == cnv$type))
})

test_that("Can write header", {
    header<- make_header(gbuild= 'hg38', genomes= genomes, is_chrom_prefixed= TRUE, cmd= '##Command foo bar', extra= c(purity= 0.5, ploidy= 3.1))
    header<- unlist(strsplit(header, split= '\n'))
    expect_equal(48, length(header))
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
    exec_snp_pileup(20, 'data/common.sample.vcf.gz', 'tmp_testthat/tmp.cvs.gz', 
        normal_bam= 'data/TCRBOA6-N-WEX.sample.bam', 
        tumour_bam= 'data/TCRBOA6-T-WEX.sample.bam', 
        mapq= 10, 
        baq= 10, 
        pseudo_snp= 165)
    expect_true(file.exists('tmp_testthat/tmp.cvs.gz'))
    expect_equal(0, system2(c('gzip', '--test', 'tmp_testthat/tmp.cvs.gz')))

    csv<- read.table('tmp_testthat/tmp.cvs.gz', sep= ',', header= TRUE)
    expect_equal(12, ncol(csv))
    expect_true(nrow(csv) > 1000)
    expect_equal(20, unique(csv$Chromosome)) 

    # Coverage at this position. Checked with:
    # samtools mpileup data/TCRBOA6-T-WEX.sample.bam -r 20:25058550 | head
    expect_equal(54, csv[csv$Position ==  25058550, "File1R"])
    expect_equal(0, csv[csv$Position ==  25058550, "File1A"])
    expect_equal(65, csv[csv$Position ==  25058550, "File2R"])
    expect_equal(0, csv[csv$Position ==  25058550, "File2A"])
})

test_that("Can execute parallel pileup", {
    mkdir()
    # snp_vcf, output, normal_bam, tumour_bam, mapq, baq, pseudo_snp, nprocs
    exec_snp_pileup_parallel('data/common.sample.vcf.gz', 'tmp_testthat/tmp.cvs.gz', 
        normal_bam= 'data/TCRBOA6-N-WEX.sample.bam', 
        tumour_bam= 'data/TCRBOA6-T-WEX.sample.bam', 
        mapq= 10, 
        baq= 10, 
        pseudo_snp= 250,
        nprocs= 3)
    expect_true(file.exists('tmp_testthat/tmp.cvs.gz'))
    csv<- fread('gzip -c -d tmp_testthat/tmp.cvs.gz')
    expect_equal(3, length(unique(csv$Chromosome)))
})

test_that("Pileup throws error if bash script throws error", {
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

test_that("Can plot coverage histogram", {
    mkdir()
    rcmat<- fread('gzip -c -d data/rcmat.txt.gz')
    plot_coverage(rcmat, rcmat[NOR.DP > 100], fname= 'tmp_testthat/hist.pdf', title= 'Depth of coverage\nmy/output/prefix')
    expect_true(file.size('tmp_testthat/hist.pdf') > 1000)
})

unlink('tmp_testthat', recursive= TRUE)
