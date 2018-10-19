#!/usr/bin/env bats
#
# The bats testing framework can be found at
# https://github.com/sstephenson/bats and it should have been installed by
# setup.sh 
#
teardown() {
    rm -rf test_out
}

@test "Can show help" {
    run ../bin/cnv_facets.R --help
    [ "$status" -eq 0 ]
    [[ "$output" == *"--tumour"* ]]
}

@test "Can use compressed pileup input" {
    run ../bin/cnv_facets.R -p data/pileup.csv.gz -o test_out/out
    [ "$status" -eq 0 ]
    run ls test_out/out.vcf.gz test_out/out.cnv.png test_out/out.spider.pdf
    [ "$status" -eq 0 ]
}

@test "Can use uncompressed pileup input" {
    run ../bin/cnv_facets.R -p data/pileup.csv -o test_out/out
    [ "$status" -eq 0 ]
    run ls test_out/out.vcf.gz test_out/out.cnv.png test_out/out.spider.pdf
    [ "$status" -eq 0 ]
    run bash -c "zcat test_out/out.vcf.gz | grep -q -F '##contig=<ID=chr1>'"
    [ "$status" -eq 0 ]
}

@test "Can use BAM input" {
    run ../bin/cnv_facets.R -t data/tumour.bam -n data/normal.bam -vcf data/snps.vcf.gz -o test_out/out
    [ "$status" -eq 0 ]
    run ls test_out/out.vcf.gz test_out/out.cnv.png test_out/out.spider.pdf
    [ "$status" -eq 0 ]
}

@test "Run real dataset" {
    run ../bin/cnv_facets.R -p data/stomach.csv.gz -o test_out/out
    [ "$status" -eq 0 ]
    run ls test_out/out.vcf.gz test_out/out.cnv.png test_out/out.spider.pdf
    [ "$status" -eq 0 ]
}
