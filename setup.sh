#!/usr/bin/env bash
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

VERSION=0.2.0

set -e
set -o pipefail

# Parse arguments
# ===============

PG=`basename "$0"`
bin_dir=${HOME}/bin

# From https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -b|--bin_dir)
        bin_dir="$2"
        shift # past argument
        shift # past value
    ;;
    -xt|--skip_test)
        skip_test=1
        shift 
    ;;
    -v|--version)
        version=1
        shift 
    ;;
    -h|--help)
        help=1
        shift
    ;;
    *)
        echo "Unknown option in $@"  
        exit 1
    shift # past argument
    ;;
esac
done

if [[ $help == 1 ]]
then
cat <<EOF
DESCRIPTION
Installer of several dependencies for cnv_facets.R. NB: Dependencies already
found are not re-installed but their version is not checked.

-b|--bin_dir     Install missing programs here. This dir should writable and on
                 your PATH. Default $bin_dir
-xt|--skip_test  Do not run the test suite at the end of the installation
-v|--version     Show version
-h|--help        Show help

USAGE EXAMPLE
bash setup.sh -b $bin_dir

Version $VERSION
EOF
exit 0
fi

if [[ $version == 1 ]]
then
    echo "$PG $VERSION"
    exit 0
fi

# End argument parsing
# ====================

function check_bin_dir(){
    # USAGE: check_bin_dir /foo/bar
    # 
    # Check `/foo/bar` is on PATH. If it is on PATH but it does not exist,
    # try to create it
    python3 -c "import os, sys
PATH= os.environ['PATH'].split(os.path.pathsep)
PATH= [os.path.abspath(x) for x in PATH]
bin= os.path.abspath('$1')
if bin not in PATH:
    sys.stderr.write('\n\033[31mError: bin directory requested by -b/--bin_dir \'$1\' is not on PATH\033[0m\n')
    sys.exit(1)
if not os.path.isdir(bin):
    os.makedirs(bin) 
"
}

function install_htslib(){
    # Download and install htslib. Compiled stuff is in `pwd`/htslib 
    pushd .
    rm -f htslib-1.8.tar.bz2
    wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2
    tar xf htslib-1.8.tar.bz2
    rm htslib-1.8.tar.bz2
    mv htslib-1.8 htslib
    cd htslib
    ./configure --prefix=`pwd`
    make -j 4
    make install
    popd 
}

# =======================

check_bin_dir ${bin_dir}

cwd=`pwd`
mkdir -p tmp

# R packages
cd ${cwd}
Rscript install/install_pkgs.R

# Samtools
found=`command -v samtools` || true
if [[ -z $found ]]
then
    cd ${cwd}/tmp
    rm -f samtools-1.8.tar.bz2
    wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2
    tar xf samtools-1.8.tar.bz2
    rm samtools-1.8.tar.bz2
    cd samtools-1.8
    ./configure --prefix=`pwd`
    make -j 4
    cp samtools ${bin_dir}/
fi
command -v samtools
samtools --version

# bcftools
found=`command -v bcftools` || true
if [[ -z $found ]]
then
    cd ${cwd}/tmp
    rm -f bcftools-1.8.tar.bz2
    wget https://github.com/samtools/bcftools/releases/download/1.8/bcftools-1.8.tar.bz2
    tar xf bcftools-1.8.tar.bz2
    rm bcftools-1.8.tar.bz2
    cd bcftools-1.8
    ./configure --prefix=`pwd`
    make -j 4
    cp bcftools ${bin_dir}/
fi
command -v bcftools
bcftools --version

# FACETS::snp-pileup
found=`command -v snp-pileup` || true
if [[ -z $found ]]
then
    cd ${cwd}/tmp
    rm -rf facets
    git clone https://github.com/mskcc/facets.git
    cd facets/inst/extcode/
    install_htslib
    ln -s `pwd`/htslib/lib/libhts.a `pwd`/htslib/lib/libhts-static.a
    g++ -std=c++11 -I `pwd`/htslib/include snp-pileup.cpp -L `pwd`/htslib/lib -lhts-static -o snp-pileup -lcurl -lz -lpthread -lcrypto -llzma -lbz2
    cp `pwd`/snp-pileup ${bin_dir}/
    rm -rf `pwd`/htslib
fi
command -v snp-pileup
snp-pileup --help

# FACETS::Rscript
cd ${cwd}
chmod a+x bin/cnv_facets.R
rsync --update bin/cnv_facets.R ${bin_dir}
cnv_facets.R --help

# VCF validator - Only for testing
if [[ "$OSTYPE" == "darwin"* ]]; 
then
    curl -L -o test/vcf_validator -s https://github.com/EBIvariation/vcf-validator/releases/download/v0.9.1/vcf_validator_macos
else
    curl -L -o test/vcf_validator -s https://github.com/EBIvariation/vcf-validator/releases/download/v0.9.1/vcf_validator_linux
fi
chmod a+x test/vcf_validator
./test/vcf_validator --help

set +x
rm -rf tmp
echo -e "\n\033[32mSetup successfully completed\033[0m"

# TEST SUITE

if [[ $skip_test != 1 ]]
then
    echo -e "\n\033[1mExecuting test suite...\033[0m\n"
    cd ${cwd}/test
    chmod a+x test_cnv_facets.py
    ./test_cnv_facets.py
fi

