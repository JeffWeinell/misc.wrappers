% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_fastStructure.R
\name{run_fastStructure}
\alias{run_fastStructure}
\title{Run run_fastStructure from SNP data in a VCF file and and plot results.

This fucntion is a wrapper that enables running fastStructure with SNP data in either a VCF file or vcfR object and coordinates in file or a matrix or data frame object.

Notes for running fastStructure:
Requires python 2
fastStructure 'cv' (cross-validation) feature is not yet implemented here.}
\usage{
run_fastStructure(
  x,
  format = "VCF",
  coords = NULL,
  samplenames = NULL,
  kmax = 10,
  save.in = NULL,
  reps = 30,
  tolerance = 1e-05,
  prior = "simple",
  full = FALSE,
  seed = NULL,
  python.path = NULL,
  fastStructure.path = NULL,
  cleanup = TRUE,
  include.out = c(".pdf", ".Qlog", ".margLlog", ".extraLog", ".Plog"),
  debug = FALSE,
  overwrite = FALSE
)
}
\arguments{
\item{x}{'vcfR' object (see package::vcfR) or a character string with path to a SNPs dataset formatted according to the 'format' argument. Currently VCF or 'fastStructure' (a type of STRUCTURE format) can be used.}

\item{format}{Character string indicating the format of the data. Currently only "VCF" or "fastStructure" allowed. Other types may be added. Ignored if x is a vcfR object.}

\item{coords}{Either a character string with path to file containing coordinates (longitude in first column, latitude in second column), or matrix object with longitude and latitude columns.}

\item{samplenames}{NULL or a character string vector with names of samples in the input data, and coords file if supplied. If NULL (the default), sample names are extracted from the SNPs datafile.}

\item{kmax}{Numerical vector with set of values to use for K. Default 40.}

\item{save.in}{Character string with path to directory where output files should be saved.}

\item{reps}{Number of repititions. Default 100.}

\item{tolerance}{Tolerance for convergence, i.e., the change in marginal likelihood required to continue.}

\item{prior}{Type of prior to use. Default "simple".}

\item{full}{Whether or not to generate output files holding variation of Q, P, and marginal likelihood, in addition to the files holding means. Default FALSE.}

\item{seed}{Value to use as a seed for reproducing results. Default NULL.}

\item{python.path}{Character string with path to python 2 with fastStructure dependencies Numpy, Scipy, Cython, GNU Scientific Library}

\item{fastStructure.path}{Character string with path to folder containing the fastStructure python executable called 'structure.py'}

\item{cleanup}{Whether or not the original fastStructure output files (*.log, *.meanQ, *meanP file for each replicate of each K) should be deleted after the data from those files are compiled and saved in three tables. Default TRUE.}

\item{include.out}{Character vector indicating which type of files should be included as output. Default is c(".pdf",".Qlog",".margLlog"). An additional file ".Plog" can be included but can be very large.}

\item{debug}{Logical indicating whether or not to print messages indicating the internal step of the function.}

\item{overwrite}{Logical indicating whether or not to allow new output files to overwrite existing ones. Default FALSE.}

\item{save.as}{Where to save the output PDF. Default is NULL.}
}
\value{
List of plots
}
\description{
Run run_fastStructure from SNP data in a VCF file and and plot results.

This fucntion is a wrapper that enables running fastStructure with SNP data in either a VCF file or vcfR object and coordinates in file or a matrix or data frame object.

Notes for running fastStructure:
Requires python 2
fastStructure 'cv' (cross-validation) feature is not yet implemented here.
}
