% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_eems.R
\name{eems_contours}
\alias{eems_contours}
\title{eems contours

Internal function used for making EEMS plots.
This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').}
\usage{
eems_contours(mcmcpath, dimns, longlat, plot_params, is_mrates)
}
\arguments{
\item{mcmcpath}{Path the directory containing EEMs output files for a particular chain}

\item{dimns}{output of read_dimns function}

\item{longlat}{Whether or not to use long lat coordinates. Default TRUE}

\item{plot_params}{List of plotting parameters}

\item{is_mrates}{Logical indicating whether or not to generated the first two EEMs plots or the third and fourth EEMs plot}
}
\value{
A list with parameters to use for generating EEMS plots.
}
\description{
eems contours

Internal function used for making EEMS plots.
This function is from the R package 'reemsplots2' by Desislava Petkova (github repository 'dipetkov/reemsplots2').
}
