% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stacks.R
\name{tsv2bam_setup}
\alias{tsv2bam_setup}
\title{batch setup tsv2bam

For each set samples (processed with ustacks, cstacks, and sstacks) and a population map, this function creates a bash job file to run tsv2bam.
The output of this function is a character string vector with paths to bash job files, which can be piped to the function 'rsbatch' to submit the jobs.}
\usage{
tsv2bam_setup(
  save.in,
  popmap = NULL,
  samples = NULL,
  samples.dir = NULL,
  t = 16,
  tsv2bam.path = NULL
)
}
\arguments{
\item{save.in}{Character string with directory where output directories and files should be saved. This should be the same as the value used for the 'save.in' argument of ustacks_setup.}

\item{popmap}{Character string to population map. Default NULL, in which case all individuals in 'save.in' are assigned to the same 'population'.}

\item{samples}{Character string vector with paths to input '.fq.gz' data files (i.e., those supplied to ustacks). Alternative to using 'samples.dir' argument. Default NULL. If 'samples' is also NULL, then data are assumed to be single rather than paired end reads.}

\item{samples.dir}{Character string with directory containing '.fq.gz' data files. Alternative to using 'samples' argument. Default NULL. If 'samples' is also NULL, then data are assumed to be single rather than paired end reads.}

\item{t}{Number of threads to use for job parallization. Default 16.}

\item{tsv2bam.path}{Character string with path to the tsv2bam executable. Default NULL, in which case the path is derived from the default setting, which can be set with config_miscwrappers('<path/to/tsv2bam>').}
}
\value{
Character string vector with the paths to bash job files.
}
\description{
batch setup tsv2bam

For each set samples (processed with ustacks, cstacks, and sstacks) and a population map, this function creates a bash job file to run tsv2bam.
The output of this function is a character string vector with paths to bash job files, which can be piped to the function 'rsbatch' to submit the jobs.
}
