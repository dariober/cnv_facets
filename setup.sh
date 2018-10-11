#!/usr/bin/env bash
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

VERSION=0.1.0

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
    -v|--version)
        version=1
        shift # past argument
    ;;
    -h|--help)
        help=1
        shift # past argument
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
Installer of several FACETS dependencies.

-b|--bin_dir  Install missing programs here. This dir should writable and on
              your PATH. Default $bin_dir
-v|--version  Show version
-h|--help     Show help

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

cwd=`pwd`
mkdir -p tmp

# R packages
cd ${cwd}
Rscript install/install_pkgs.R

# FACETS::snp-pileup
found=`command -v snp-pileup` || true
if [[ -z $found ]]
then
    cd ${cwd}/tmp
    rm -rf facets
    git clone https://github.com/mskcc/facets.git
    cd facets/inst/extcode/
    install_htslib
    g++ -std=c++11 -I`pwd`/htslib/include snp-pileup.cpp \
        -L`pwd`/htslib/lib -lhts -Wl,-rpath=`pwd`/htslib/lib -o snp-pileup
    ln -sf `pwd`/snp-pileup ${bin_dir}/
fi
command -v snp-pileup
snp-pileup --help

# FACETS::Rscript
cd ${cwd}
chmod a+x bin/cnv_facets.R
rsync --update bin/cnv_facets.R ${bin_dir}
cnv_facets.R --help

set +x
echo -e "\n\033[32mSetup successfully completed\033[0m\n"
