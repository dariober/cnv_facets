#!/usr/bin/env python
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
import subprocess as sp
import sys
import gzip
import re
import filecmp
import gzip
from collections import OrderedDict

def vcf_to_list(vcf_file):
    vcf= []
    with gzip.open(vcf_file) as gz:
        for line in gz:
            vcf.append(line.decode())
    return vcf

def vcf_validator(vcf_file):
    """Run vcf_validator and scan the report file. Return empty string if vcf
    file is valid otherwise print the content of the report
    """
    p = sp.Popen("gzip -c -d %s | ./vcf_validator --report text" % vcf_file, shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
    stdout, stderr = p.communicate()
    log= stderr.decode().split('\n')
    rfile= [x for x in log if x.startswith('[info] Text report written to : ')]
    assert len(rfile) == 1
    rfile= rfile[0].replace('[info] Text report written to : ', '')
    with open(rfile) as f:
        report= f.readlines()
    os.remove(rfile)
    if len(report) == 1 and 'the input file is valid' in report[0]:
        return ''
    return ''.join(''.join(report))

class cnv_facets(unittest.TestCase):

    def setUp(self):
        sys.stderr.write('\n' + self.id().split('.')[-1] + ' ') # Print test name
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')

    def tearDown(self):
        if os.path.exists('test_out'):
            shutil.rmtree('test_out')

    def testShowHelp(self):
        p = sp.Popen("../bin/cnv_facets.R --help", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue('--tumour', stderr.decode())
        
    def testBamInput(self):
        p = sp.Popen("../bin/cnv_facets.R -d 1 8000 -t data/TCRBOA6-T-WEX.sample.bam -n data/TCRBOA6-N-WEX.sample.bam -vcf data/common.sample.vcf.gz -o test_out/out", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.exists('test_out/out.vcf.gz'))
        self.assertTrue(os.path.exists('test_out/out.cnv.png'))
        self.assertTrue(os.path.exists('test_out/out.cov.pdf'))
        self.assertTrue(os.path.exists('test_out/out.spider.pdf'))
        self.assertTrue(os.path.exists('test_out/out.csv.gz'))
        self.assertEqual('', vcf_validator('test_out/out.vcf.gz'))

    def testPileupInput(self):
        p = sp.Popen("../bin/cnv_facets.R --pileup data/stomach.csv.gz -o test_out/out", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.exists('test_out/out.vcf.gz'))
        self.assertTrue(os.path.exists('test_out/out.cnv.png'))
        self.assertTrue(os.path.exists('test_out/out.cov.pdf'))
        self.assertTrue(os.path.exists('test_out/out.spider.pdf'))
        self.assertEqual('', vcf_validator('test_out/out.vcf.gz'))

    def testFailOnSnpPileup(self):
        p = sp.Popen("../bin/cnv_facets.R -t data/INVALID.bam -n data/TCRBOA6-N-WEX.sample.bam -vcf data/common.sample.vcf.gz -o test_out/out", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertTrue(p.returncode != 0)
        # Check we exited immediatly after failing the first snp-pileup
        self.assertEqual(1, stderr.count('samtools view: failed to open'))

    def testSingleEndBam(self):
        # Prepare single-end files
        p= sp.check_call(r"""
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
        """, shell= True)
        
        p = sp.Popen("../bin/cnv_facets.R -t test_out/TCRBOA6-T-WEX.se.bam -n test_out/TCRBOA6-N-WEX.se.bam -vcf data/common.sample.vcf.gz -o test_out/out", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.exists('test_out/out.csv.gz'))
        self.assertEqual('', vcf_validator('test_out/out.vcf.gz'))
        self.assertTrue(os.path.getsize("test_out/out.csv.gz") > 60000 and 
                os.path.getsize("test_out/out.csv.gz") < 70000)

        p = sp.Popen("../bin/cnv_facets.R -t test_out/TCRBOA6-T-WEX.se.bam --nbhd-snp 10 -n test_out/TCRBOA6-N-WEX.se.bam -vcf data/common.sample.vcf.gz -o test_out/out", shell=True, stdout= sp.PIPE, stderr= sp.PIPE)
        stdout, stderr = p.communicate()
        self.assertEqual(0, p.returncode)
        self.assertTrue(os.path.exists('test_out/out.csv.gz'))
        self.assertEqual('', vcf_validator('test_out/out.vcf.gz'))
        self.assertTrue(os.path.getsize("test_out/out.csv.gz") > 300000 and 
                os.path.getsize("test_out/out.csv.gz") < 400000)

if __name__ == '__main__':
    unittest.main()
