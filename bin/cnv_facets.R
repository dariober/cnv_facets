#!/usr/bin/env Rscript
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

suppressMessages(library(argparse))
suppressMessages(library(facets))
suppressMessages(library(data.table))

# -----------------------------------------------------------------------------

VERSION= sprintf('0.7.1; facets %s', packageVersion('facets'))

docstring<- sprintf('DESCRIPTION \\n\\
Analyse input pileup file for CNV using facets (https://github.com/mskcc/facets). \\n\\
See facets docs for details of input and output. \\n\\
\\n\\
The estimates of purity, ploidy and sample information from facets are in the \\n\\
header. The REF allele in output is always N. \\n\\
\\n\\
Version %s', VERSION)

parser<- ArgumentParser(description= docstring, formatter_class= 'argparse.RawTextHelpFormatter')

parser$add_argument('--pileup', '-p', help= 'Pileup for matched normal (first sample) and tumour (second sample). \\n\\
This is the output of the snp-pileup program accompanying with facets.', required= TRUE)

def<- '-'
parser$add_argument('--out', '-o', help= sprintf('Output file in VCF format. If ending in .gz it will be compressed and indexed. \\n\\
If - write to stdout. Default %s', def), required= FALSE, default= def)

parser$add_argument('--plot_cnv', '-plot_cnv', help= 'Output file for CNV profile. It will be in PNG format', required= FALSE)
parser$add_argument('--plot_spider', '-plot_spider', help= 'Output for the diagnostic spider plot. It will be in PDF format', required= FALSE)

def<- 25
parser$add_argument('--ndepth', '-d', help= sprintf('Minimum depth in normal sample for a position to be considered. Default %s', def), type= 'integer', default= def)

def<- 150 
parser$add_argument('--cval', '-cval', help= sprintf('Critical value for segmentation. Default %s', def), type= 'double', default= def)

def<- 250 
parser$add_argument('--snp_nbhd', '-snp', help= sprintf('Minimum spacing between SNP before they are sampled. \\n\\
Sufficient spacing is required to reduce serial correlation. \\n\\
This option passed to prepProcSample. Default %s', def), type= 'integer', default= def)

parser$add_argument('--fai', '-f', help= sprintf('Tab separated file listing reference contigs (1st column) and their size (2nd column). \\n\\
The genome fasta index .fai is suitable. This option is used to write contigs \\n\\
in the header, not required but recommended.'), type= 'character', required= FALSE)

def<- 'hg38'
parser$add_argument('--gbuild', '-g', help= sprintf('String indicating the reference genome build. Chromosomes used for analysis are \\n\\
1-22, X for humans and 1-19 for mouse. Default %s.', def), 
    type= 'character', required= FALSE, default= def, 
    choices= c('hg18', 'hg19', 'hg38', 'mm9', 'mm10'))

def_rnd<- 'The name of the input file'
parser$add_argument('--rnd_seed', '-s', help= sprintf('Seed for random number generator. Default: %s', def_rnd), type= 'character', default= def_rnd)

# NB: argparse v1.1.1+ required for -v option to work.
parser$add_argument("-v", "--version", action= 'version', version= VERSION)

xargs<- parser$parse_args()

# -----------------------------------------------------------------------------
if(xargs$rnd_seed == def_rnd){
    seed<- sum(utf8ToInt(xargs$pileup))
} else if( ! is.na(suppressWarnings(as.numeric(xargs$rnd_seed)))){
    seed<- as.numeric(xargs$rnd_seed)
} else {
    seed<- sum(utf8ToInt(xargs$rnd_seed))
}
set.seed(seed)

readSnpMatrix2<- function(pileup){
    rcmat<- fread(sprintf('gunzip -c %s', pileup), select= c('Chromosome', 'Position', 'File1R', 'File1A', 'File2R', 'File2A'))
    setnames(rcmat, c('File1R', 'File1A', 'File2R', 'File2A'),
                    c('NOR.RD', 'NOR.DP', 'TUM.RD', 'TUM.DP'))
    rcmat[, NOR.DP := NOR.DP + NOR.RD]
    rcmat[, TUM.DP := TUM.DP + TUM.RD]
    rcmat[, Chromosome := sub("chr", "", Chromosome, fixed= TRUE)]
    setcolorder(rcmat, c('Chromosome', 'Position', 'NOR.DP', 'NOR.RD', 'TUM.DP', 'TUM.RD'))
    return(rcmat)
}

facetsRecordToVcf<- function(x){
    # Convert the annotated facets record to a VCF record.
    vcf<- vector(length= 8)
    vcf[1]<- x$chrom
    vcf[2]<- x$start+1 # TODO: Check +1 is correct
    vcf[3]<- x$seg
    vcf[4]<- 'N'
    vcf[5]<- '<CNV>'
    vcf[6]<- '.'
    vcf[7]<- ifelse(x$type == "NEUTR", 'neutral', 'PASS')

    # INFO field: Keep consistent with header
    vcf[8]<- paste0(
        'SVTYPE=', x$type,
        ';SVLEN=', x$end-x$start,
        ';END=', x$end,
        ';NUM_MARK=', x$num.mark,
        ';NHET=', x$nhet,
        ';CNLR_MEDIAN=', ifelse(is.na(x$cnlr.median), '.', round(x$cnlr.median, 3)), 
        ';MAF_R=', ifelse(is.na(x$mafR), '.', round(x$mafR, 3)), 
        ';SEGCLUST=', x$segclust, 
        ';CNLR_MEDIAN_CLUST=', ifelse(is.na(x$cnlr.median.clust), '.', round(x$cnlr.median.clust, 3)), 
        ';MAF_R_CLUST=', ifelse(is.na(x$mafR.clust), '.', round(x$mafR.clust, 3)), 
        ';CF_EM=', ifelse(is.na(x$cf.em), '.', round(x$cf.em, 3)),
        ';TCN_EM=', ifelse(is.na(x$tcn.em), '.', x$tcn.em),
        ';LCN_EM=', ifelse(is.na(x$lcn.em), '.', x$lcn.em))
   return(vcf)
}

getScriptName<- function(){
    opt<- grep('^--file=', commandArgs(trailingOnly = FALSE), value= TRUE)
    name<- basename(sub('^--file=', '', opt))
    return(name)
}

write(sprintf('Loading file %s...', xargs$pileup), stderr())
rcmat<- readSnpMatrix2(xargs$pileup)

write(sprintf('Preprocessing sample...'), stderr())
xx<- preProcSample(rcmat, ndepth= xargs$ndepth, gbuild= xargs$gbuild, snp.nbhd= xargs$snp_nbhd, het.thresh=0.25, cval=25, deltaCN=0, unmatched= FALSE, ndepthmax= 1000)
rm(rcmat)
x_ <- gc(verbose= FALSE)

write(sprintf('Processing sample...'), stderr())
oo<- procSample(xx, cval= xargs$cval, min.nhet=15, dipLogR=NULL)

write(sprintf('Fitting model...'), stderr())
fit<- emcncf(oo, unif= FALSE, min.nhet= 15, maxiter= 10, eps=1e-3)

write(sprintf('Writing output'), stderr())
out<- data.table(fit$cncf)
stopifnot(all(!grepl('chr', out$chrom)))
out[, chrom := paste0('chr', chrom)]
out[, chrom := sub('chr23', 'chrX', chrom)] # This is because of bug https://github.com/mskcc/facets/issues/60
setcolorder(out, c('chrom', 'start', 'end', 'seg', 'num.mark', 'nhet', 'cnlr.median', 'mafR', 'segclust', 'cnlr.median.clust', 'mafR.clust', 'cf.em', 'tcn.em', 'lcn.em'))

# Classify CNV. See also https://github.com/mskcc/facets/issues/62
out[, type := NA]
out[, type := ifelse((tcn.em == 2 & (lcn.em == 1 | is.na(lcn.em))), 'NEUTR', type)]
out[, type := ifelse(is.na(type) & tcn.em == 0, 'DEL', type)]
out[, type := ifelse(is.na(type) & tcn.em > 2 & (lcn.em > 0 | is.na(lcn.em)), 'DUP', type)]
out[, type := ifelse(is.na(type) & tcn.em == 1, 'HEMIZYG', type)]
out[, type := ifelse(is.na(type) & tcn.em == 2 & lcn.em == 0, 'LOH', type)]
out[, type := ifelse(is.na(type) & tcn.em > 2 & lcn.em == 0, 'DUP-LOH', type)]
out<- out[order(chrom, start)]
stopifnot(all(!is.na(out$type))) # Everything has been classified

if(xargs$out == '-'){
    xargs$out<- stdout()
}

write('##fileformat=VCFv4.2', xargs$out, append= FALSE)
write('##FILTER=<ID=PASS,Description="All filters passed">', xargs$out, append= TRUE)
write(sprintf('##%sCommand=%s; Version=%s; Date=%s', getScriptName(), paste(commandArgs(), collapse= ' '), VERSION, Sys.time()), xargs$out, append= TRUE)
if( is.null(xargs$fai) ){
    for(chrom in unique(out$chrom)){
        write(sprintf('##contig=<ID=%s>', chrom), xargs$out, append= TRUE)
    }
} else{
    fai<- fread(xargs$fai)
    for(i in 1:nrow(fai)){
        write(sprintf('##contig=<ID=%s,length=%s>', fai[i]$V1, fai[i]$V2), xargs$out, append= TRUE)
    }
}
write('##FILTER=<ID=PASS,Description="All filters passed">', xargs$out, append= TRUE)
write('##FILTER=<ID=neutral,Description="Copy number neutral">', xargs$out, append= TRUE)
write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">', xargs$out, append= TRUE)
write('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">', xargs$out, append= TRUE)
write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">', xargs$out, append= TRUE)
write('##INFO=<ID=NUM_MARK,Number=1,Type=Integer,Description="Number of SNPs in the segment">', xargs$out, append= TRUE)
write('##INFO=<ID=NHET,Number=1,Type=Integer,Description="Number of SNPs that are deemed heterozygous">', xargs$out, append= TRUE)
write('##INFO=<ID=CNLR_MEDIAN,Number=1,Type=Float,Description="Median log-ratio (logR) of the segment. logR is defined by the log-ratio of total read depth in the tumor versus that in the normal">', xargs$out, append= TRUE)
write('##INFO=<ID=CNLR_MEDIAN_CLUST,Number=1,Type=Float,Description="Median log-ratio (logR) of the segment cluster. logR is defined by the log-ratio of total read depth in the tumor versus that in the normal">', xargs$out, append= TRUE)
write('##INFO=<ID=MAF_R,Number=1,Type=Float,Description="Log-odds-ratio (logOR) summary for the segment. logOR is defined by the log-odds ratio of the variant allele count in the tumor versus in the normal">', xargs$out, append= TRUE)
write('##INFO=<ID=MAF_R_CLUST,Number=1,Type=Float,Description="Log-odds-ratio (logOR) summary for the segment cluster. logOR is defined by the log-odds ratio of the variant allele count in the tumor versus that in the normal">', xargs$out, append= TRUE)
write('##INFO=<ID=SEGCLUST,Number=1,Type=Integer,Description="Segment cluster to which the segment belongs">', xargs$out, append= TRUE)
write('##INFO=<ID=CF_EM,Number=1,Type=Float,Description="Cellular fraction, fraction of DNA associated with the aberrant genotype. Set to 1 for normal diploid">', xargs$out, append= TRUE)
write('##INFO=<ID=TCN_EM,Number=1,Type=Integer,Description="Total copy number. 2 for normal diploid">', xargs$out, append= TRUE)
write('##INFO=<ID=LCN_EM,Number=1,Type=Integer,Description="Lesser (minor) copy number. 1 for normal diploid">', xargs$out, append= TRUE)

write(paste0('##purity=', fit$purity), xargs$out, append= TRUE)
write(paste0('##ploidy=', fit$ploidy), xargs$out, append= TRUE)
write(paste0('##dipLogR=', fit$dipLogR), xargs$out, append= TRUE)
write(paste0('##emflags="', fit$emflags, '"'), xargs$out, append= TRUE)

write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO', xargs$out, append= TRUE)
for(i in 1:nrow(out)){
    write(paste(facetsRecordToVcf(out[i]), collapse= '\t'), xargs$out, append= TRUE)
}

if(grepl('\\.gz$', xargs$out)){
    tmp<- tempfile(pattern= basename(xargs$out), tmpdir= dirname(xargs$out))
    exit_code<- system(sprintf('bgzip -f -c %s > %s', xargs$out, tmp))
    stopifnot(exit_code == 0)
    rn<- file.rename(tmp, xargs$out)
    exit_code<- system(sprintf('tabix -f %s', xargs$out))
    stopifnot(exit_code == 0)
}

if( ! is.null(xargs$plot_cnv) ){
    write(sprintf('Plotting genome...'), stderr())
    png(xargs$plot_cnv, units="px", width=1600, height=1600, res=300)
    sname<- sub('.csv.gz$|.csv$', '', basename(xargs$pileup))
    sname<- sprintf('%s; ploidy= %.2f; purity= %.2f', sname, fit$ploidy, fit$purity)
    plotSample(x=oo, emfit=fit, sname= sname)
    x_ <- dev.off()
}

if( ! is.null(xargs$plot_spider) ){
    write(sprintf('Plotting spider...'), stderr())
    pdf(xargs$plot_spider, width= 16/2.54, height= 14/2.54)
    par(las= 1, mar= c(3, 3, 1, 1), mgp= c(1.5, 0.5, 0), tcl= -0.3)
    logRlogORspider(oo$out, oo$dipLogR)
    x_ <- dev.off()
}

si<- capture.output(sessionInfo())
write(si, stderr())
