% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stacks.R
\name{summarize_VCFs}
\alias{summarize_VCFs}
\title{summarize_VCFs

Same as 'summarize_stacks', except the full paths to input VCFs and an output directory are supplied as arguments.}
\usage{
summarize_VCFs(
  VCF.paths,
  save.as,
  popmap.path = NULL,
  include.out = c("total", "persite")
)
}
\arguments{
\item{VCF.paths}{Character string vector to one or more input VCF files.}

\item{save.as}{Path to location where output stats table should be saved. If NULL, the the table is not saved.}

\item{popmap.path}{Character string to population map. If NULL, all individuals in are assigned to the same 'population'.}

\item{include.out}{Character string vector indicating which output files to generate.}
}
\value{
data frame(s) with basic stats
}
\description{
summarize_VCFs

Same as 'summarize_stacks', except the full paths to input VCFs and an output directory are supplied as arguments.
}
