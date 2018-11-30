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

repos<- 'https://cran.r-project.org'

lib<- NULL
for(x in .libPaths()){
    if(file.access(x, mode= 2) == 0){
        lib<- x
        break
    }
}

if(is.null(lib)){
    stop(sprintf('Cannot find a writable directory to install R packages.'))
}

tryCatch({
        suppressMessages(library(devtools))
    }, error = function(err) {
        install.packages('devtools', repos= repos, lib= lib)
    }
)

tryCatch({
        suppressMessages(library(Rsamtools))
    }, error = function(err) {
        source("https://bioconductor.org/biocLite.R")
        biocLite("Rsamtools", ask= FALSE, suppressUpdates= TRUE, lib= lib)
    }
)

tryCatch({
        suppressMessages(library(data.table))
    }, error = function(err) {
        install.packages('data.table', repos= repos, lib)
    }
)

tryCatch({
        suppressMessages(library(argparse))
    }, error = function(err) {
        devtools::install_github("trevorld/argparse", ref= "v1.1.1", repos= repos, lib= lib)
    }
)

tryCatch({
        suppressMessages(library(facets))
    }, error = function(err) {
        devtools::install_github("mskcc/facets", ref = "434b5ce", lib= lib)
    }
)

tryCatch({
        suppressMessages(library(testthat))
    }, error = function(err) {
        install.packages('testthat', repos= repos, lib)
    }
)
