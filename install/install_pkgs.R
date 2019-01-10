#!/usr/bin/env Rscript
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

# Download missing packages from this CRAN mirror
REPOS<- 'https://cran.r-project.org'

# Install missing R packages in this directory. If NULL, attempt to find or
# create a suitable one. If not NULL, this directory must exist and be on the
# path of R libraries
LIB<- NULL

FACETS_REF<- '434b5ce'
# -----------------------------------------------------------------------------

if(is.null(LIB)){
    # Get the first writable directory, if any
    for(x in .libPaths()){
        if(file.access(x, mode= 2) == 0){
            LIB<- x
            break
        }
    }
    if(is.null(LIB)){
        # Try to create a writable directory using the same strategy as
        # install.packages
        LIB <- unlist(strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep))[1L]
        if (!file.exists(LIB)) {
            if (!dir.create(LIB, recursive = TRUE)){
                stop(gettextf("unable to create %s", sQuote(LIB)), domain = NA)
            }
            .libPaths(c(LIB, .libPaths()))
        }
    }
}
if(is.null(LIB)){
    stop(sprintf('Cannot find a writable directory to install R packages.'))
} else {
    write(sprintf('Installing any missing R package to library: "%s"', LIB), stderr())
}

# -----------------------------------------------------------------------------

tryCatch({
        suppressMessages(library(devtools))
    }, error = function(err) {
        install.packages('devtools', repos= REPOS, lib= LIB, Ncpus= 6)
    }
)

tryCatch({
        suppressMessages(library(Rsamtools))
    }, error = function(err) {
        source("https://bioconductor.org/biocLite.R")
        biocLite("Rsamtools", ask= FALSE, suppressUpdates= TRUE, lib= LIB, Ncpus= 6)
    }
)

tryCatch({
        suppressMessages(library(data.table))
    }, error = function(err) {
        install.packages('data.table', repos= REPOS, lib= LIB, Ncpus= 6)
    }
)

tryCatch({
        suppressMessages(library(ggplot2))
    }, error = function(err) {
        install.packages('ggplot2', repos= REPOS, lib= LIB, Ncpus= 6)
    }
)

tryCatch({
        suppressMessages(library(gridExtra))
    }, error = function(err) {
        install.packages('gridExtra', repos= REPOS, lib= LIB, Ncpus= 6)
    }
)

tryCatch({
        suppressMessages(library(argparse))
    }, error = function(err) {
        devtools::install_github("trevorld/argparse", ref= "v1.1.1", repos= REPOS, lib= LIB)
    }
)

tryCatch({
        suppressMessages(library(facets))
    }, error = function(err) {
        devtools::install_github("mskcc/facets", ref = FACETS_REF, lib= LIB)
    }
)

tryCatch({
        suppressMessages(library(testthat))
    }, error = function(err) {
        install.packages('testthat', repos= REPOS, lib= LIB, Ncpus= 6)
    }
)

tryCatch({
        suppressMessages(library(covr))
    }, error = function(err) {
        install.packages('covr', repos= REPOS, lib= LIB, Ncpus= 6)
    }
)
