#!/usr/bin/env python3
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

import unittest
import shutil
import os
import subprocess
import sys
import gzip
import re
import filecmp
from collections import OrderedDict


class shell:
    def __init__(self, cmd, strict=True, timeout=None):
        print(cmd)
        cmd = f"set -e; set -u; set -o pipefail\n{cmd}"
        p = subprocess.Popen(
            cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            executable="/bin/bash",
        )
        try:
            stdout, stderr = p.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            p.kill()
            sys.stderr.write(f"Error: Timeout after {timeout} seconds\n")
            stdout, stderr = p.communicate()
        self.returncode = p.returncode
        self.stdout = stdout.decode()
        self.stderr = stderr.decode()
        self.cmd = cmd
        if strict and self.returncode != 0:
            raise subprocess.SubprocessError(
                f"\nSTDOUT:\n{self.stdout}\nSTDERR:\n{self.stderr}\nEXIT CODE: {self.returncode}"
            )


def vcf_to_list(vcf_file):
    vcf = []
    with gzip.open(vcf_file) as gz:
        for line in gz:
            vcf.append(line.decode().strip().split("\t"))
    return vcf


def vcf_validator(vcf_file):
    """Run vcf_validator and scan the report file. Return empty string if vcf
    file is valid otherwise print the content of the report
    """
    p = subprocess.Popen(
        "export LC_ALL=C; gzip -c -d %s | ./vcf_validator --report text" % vcf_file,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, stderr = p.communicate()
    log = stderr.decode().split("\n")
    rfile = [x for x in log if x.startswith("[info] Text report written to : ")]
    assert len(rfile) == 1
    rfile = rfile[0].replace("[info] Text report written to : ", "")
    with open(rfile) as f:
        report = f.readlines()
    os.remove(rfile)
    if len(report) == 1 and "the input file is valid" in report[0]:
        return ""
    return "".join("".join(report))


class cnv_facets(unittest.TestCase):
    def setUp(self):
        sys.stderr.write("\n" + self.id().split(".")[-1] + "\n")  # Print test name
        if os.path.exists("test_out"):
            shutil.rmtree("test_out")

    def tearDown(self):
        if os.path.exists("test_out"):
            shutil.rmtree("test_out")

    def testShowHelp(self):
        p = shell("../bin/cnv_facets.R --help")
        self.assertTrue("--tumour", p.stderr)

    def testBamInput(self):
        shell(
            "../bin/cnv_facets.R -d 1 8000 -t data/TCRBOA6-T-WEX.sample.bam -n data/TCRBOA6-N-WEX.sample.bam -vcf data/common.sample.vcf.gz -o test_out/out"
        )
        self.assertTrue(os.path.exists("test_out/out.vcf.gz"))
        self.assertTrue(os.path.exists("test_out/out.cnv.png"))
        self.assertTrue(os.path.exists("test_out/out.cov.pdf"))
        self.assertTrue(os.path.exists("test_out/out.spider.pdf"))
        self.assertTrue(os.path.exists("test_out/out.csv.gz"))
        self.assertEqual("", vcf_validator("test_out/out.vcf.gz"))

    def testDoNotPlot(self):
        shell(
            """../bin/cnv_facets.R --no-cov-plot -t data/TCRBOA6-T-WEX.sample.bam -n data/TCRBOA6-N-WEX.sample.bam -vcf data/common.sample.vcf.gz -o test_out/out"""
        )
        self.assertTrue(os.path.exists("test_out/out.vcf.gz"))
        self.assertTrue(os.path.exists("test_out/out.cnv.png"))
        self.assertTrue(os.path.exists("test_out/out.spider.pdf"))
        self.assertTrue(not os.path.exists("test_out/out.cov.pdf"))

    def testBamInputNotProperlyPaired(self):
        # Prepare a bam files with 'properly paired' flag removed
        shell(
            r"""
samtools view -h data/TCRBOA6-N-WEX.sample.bam \
| awk -v FS='\t' -v OFS='\t' '$1 ~ "^@" || $2 ~ "99|83|163|147" {if($1 ~ "^@"){print $0} else {$2=$2-2; print $0}}' \
| samtools view -b - > tmp.n.bam &&
samtools index tmp.n.bam &&

samtools view -h data/TCRBOA6-T-WEX.sample.bam \
| awk -v FS='\t' -v OFS='\t' '$1 ~ "^@" || $2 ~ "99|83|163|147" {if($1 ~ "^@"){print $0} else {$2=$2-2; print $0}}' \
| samtools view -b - > tmp.t.bam &&
samtools index tmp.t.bam 
                """
        )

        p = shell(
            r"""
../bin/cnv_facets.R -d 1 8000 -t tmp.t.bam -n tmp.n.bam -vcf data/common.sample.vcf.gz -o test_out/out
                """,
            strict=False,
        )
        self.assertEqual(1, p.returncode)
        self.assertTrue(os.path.exists("test_out/out.csv.gz"))
        os.remove("test_out/out.csv.gz")

        shell(
            r"""
../bin/cnv_facets.R --snp-count-orphans -d 1 8000 -t tmp.t.bam -n tmp.n.bam -vcf data/common.sample.vcf.gz -o test_out/out
                """
        )
        self.assertTrue(os.path.exists("test_out/out.csv.gz"))
        self.assertTrue(os.path.getsize("test_out/out.csv.gz") > 60000)
        self.assertTrue(os.path.getsize("test_out/out.vcf.gz") > 1000)

        for x in ["tmp.n.bam", "tmp.n.bam.bai", "tmp.t.bam", "tmp.t.bam.bai"]:
            os.remove(x)

    def testPileupInput(self):
        shell("../bin/cnv_facets.R --pileup data/stomach.csv.gz -o test_out/out")
        self.assertTrue(os.path.exists("test_out/out.vcf.gz"))
        self.assertTrue(os.path.exists("test_out/out.cnv.png"))
        self.assertTrue(os.path.exists("test_out/out.cov.pdf"))
        self.assertTrue(os.path.exists("test_out/out.spider.pdf"))
        self.assertEqual("", vcf_validator("test_out/out.vcf.gz"))

    def testFailOnSnpPileup(self):
        p = shell(
            "../bin/cnv_facets.R -t data/INVALID.bam -n data/TCRBOA6-N-WEX.sample.bam -vcf data/common.sample.vcf.gz -o test_out/out",
            strict=False,
        )
        self.assertTrue(p.returncode != 0)
        # Check we exited immediatly after failing the first snp-pileup
        self.assertEqual(1, p.stderr.count("samtools view: failed to open"))

    def testTargetPanel(self):
        shell(
            "../bin/cnv_facets.R -T data/stomach_targets_chr.bed -p data/stomach_chr.csv.gz -o test_out/out"
        )
        self.assertEqual("", vcf_validator("test_out/out.vcf.gz"))
        vcf = vcf_to_list("test_out/out.vcf.gz")
        vcf = [line for line in vcf if not line[0].startswith("#")]
        chroms = set([line[0] for line in vcf])
        self.assertTrue("chr1" in chroms)
        self.assertTrue("chr2" in chroms)
        self.assertTrue("chr3" in chroms)
        self.assertTrue("chr4" in chroms)
        self.assertTrue("chrX" in chroms)
        self.assertTrue("chr11" not in chroms)
        self.assertTrue("chr12" not in chroms)
        self.assertTrue("chr13" not in chroms)

        chr1 = [line for line in vcf if line[0] == "chr1"]
        self.assertTrue(int(chr1[0][1]) > 1000000)
        self.assertTrue(int(chr1[len(chr1) - 1][1]) < 100000000)

        # Without chr prefix
        shell(
            "../bin/cnv_facets.R -T data/stomach_targets.bed -p data/stomach.csv.gz -o test_out/out"
        )
        self.assertEqual("", vcf_validator("test_out/out.vcf.gz"))

        vcf = vcf_to_list("test_out/out.vcf.gz")
        vcf = [line for line in vcf if not line[0].startswith("#")]
        chroms = set([line[0] for line in vcf])
        self.assertTrue("1" in chroms)
        self.assertTrue("13" not in chroms)

        chr1 = [line for line in vcf if line[0] == "1"]
        self.assertTrue(int(chr1[0][1]) > 1000000)
        self.assertTrue(int(chr1[len(chr1) - 1][1]) < 100000000)

    def testSingleEndBam(self):
        # Prepare single-end files
        shell(
            r"""
        mkdir test_out
        for sample in TCRBOA6-N TCRBOA6-T
        do
            samtools view -h data/${sample}-WEX.sample.bam \
            | awk -v OFS='\t' -v FS='\t' '{if($1 ~ "^@"){print $0} 
                if($2 == 99) {
                    $2 = 0; print $0
                } else if($2 == 163) {
                    $2 = 16; print $0}
                }' \
            | samtools view -b - > test_out/${sample}-WEX.se.bam
             samtools index test_out/${sample}-WEX.se.bam
        done
        """
        )

        shell(
            "../bin/cnv_facets.R -t test_out/TCRBOA6-T-WEX.se.bam -n test_out/TCRBOA6-N-WEX.se.bam -vcf data/common.sample.vcf.gz -o test_out/out"
        )
        self.assertTrue(os.path.exists("test_out/out.csv.gz"))
        self.assertEqual("", vcf_validator("test_out/out.vcf.gz"))
        self.assertTrue(
            os.path.getsize("test_out/out.csv.gz") > 60000
            and os.path.getsize("test_out/out.csv.gz") < 70000
        )

        shell(
            "../bin/cnv_facets.R -t test_out/TCRBOA6-T-WEX.se.bam --nbhd-snp 10 -n test_out/TCRBOA6-N-WEX.se.bam -vcf data/common.sample.vcf.gz -o test_out/out"
        )
        self.assertTrue(os.path.exists("test_out/out.csv.gz"))
        self.assertEqual("", vcf_validator("test_out/out.vcf.gz"))
        self.assertTrue(
            os.path.getsize("test_out/out.csv.gz") > 300000
            and os.path.getsize("test_out/out.csv.gz") < 400000
        )


if __name__ == "__main__":
    unittest.main()
