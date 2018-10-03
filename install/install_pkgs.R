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


tryCatch({
        library(Rsamtools)
    }, error = function(err) {
        source("https://bioconductor.org/biocLite.R")
        biocLite("Rsamtools")
    }
)

tryCatch({
        library(data.table)
    }, error = function(err) {
        install.packages('data.table', repos='http://cran.r-project.org')
    }
)

tryCatch({
        library(argparse)
    }, error = function(err) {
        install.packages('argparse', repos='http://cran.r-project.org')
    }
)

tryCatch({
        # devtools only necessary to install other packages.
        library(devtools)
    }, error = function(err) {
        install.packages('devtools', repos='http://cran.r-project.org')
    }
)

tryCatch({
        library(facets)
    }, error = function(err) {
        devtools::install_github("mskcc/facets")
    }
)
