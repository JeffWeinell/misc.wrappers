% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exepaths_miscwrappers.R
\name{exepaths_miscwrappers}
\alias{exepaths_miscwrappers}
\title{exepaths_miscwrappers function

Function to return the paths of executables that were previously defined using the function congig_miscwrappers.}
\usage{
exepaths_miscwrappers(programs = NULL)
}
\arguments{
\item{programs}{Character string or vector with names of programs to return path settings for, or NULL, in which case all defined path settings are returned.}
}
\value{
Two column data.frame with program name and corresponding full path to the program executable file
}
\description{
exepaths_miscwrappers function

Function to return the paths of executables that were previously defined using the function congig_miscwrappers.
}
