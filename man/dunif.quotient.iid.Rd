% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_model.R
\name{dunif.quotient.iid}
\alias{dunif.quotient.iid}
\title{Probability density function for the quotient of two identically independent uniform distributions.

Probability density for z in Z ~ X1/X2, where X1 & X2 ~ U(a,b), a & b in Reals, and a < b.
The function dunif.quotient.a0 does not require X1 and X2 to be identical, but has the requirement that the lower limit of both distributions equal 0.}
\usage{
dunif.quotient.iid(z, a, b)
}
\arguments{
\item{z}{quantile}

\item{a}{Number indicating the lower limit of X1 and X2}

\item{b}{Number indicating the lower upper of X1 and X2}
}
\value{
Probability density of z in Z
}
\description{
Probability density function for the quotient of two identically independent uniform distributions.

Probability density for z in Z ~ X1/X2, where X1 & X2 ~ U(a,b), a & b in Reals, and a < b.
The function dunif.quotient.a0 does not require X1 and X2 to be identical, but has the requirement that the lower limit of both distributions equal 0.
}
