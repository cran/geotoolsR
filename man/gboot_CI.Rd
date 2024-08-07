% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gboot_CI.R
\name{gboot_CI}
\alias{gboot_CI}
\title{Bootstrap Confidence Interval}
\usage{
gboot_CI(x,alpha=0.05,digits=3)
}
\arguments{
\item{x}{object generate by functions \code{\link[geotoolsR]{gboot_block}},
\code{\link[geotoolsR]{gboot_cloud}}, \code{\link[geotoolsR]{gboot_cross}},
 \code{\link[geotoolsR]{gboot_solow}}, \code{\link[geotoolsR]{gboot_variogram}}}

\item{alpha}{significance level (Default=0.05).}

\item{digits}{number of decimal places.}
}
\value{
Invalid arguments will return an error message.
}
\description{
Provides a two-sided bootstrap confidence interval.
}
\details{
Examples of this function can be found in \code{\link[geotoolsR]{gboot_block}},
 \code{\link[geotoolsR]{gboot_cloud}}, \code{\link[geotoolsR]{gboot_cross}},
  \code{\link[geotoolsR]{gboot_solow}}, \code{\link[geotoolsR]{gboot_variogram}}
}
\author{
Diogo Francisco Rossoni \email{dfrossoni@uem.br}

Vinicius Basseto Felix \email{felix_prot@hotmail.com}
}
\keyword{Bootstrap}
\keyword{CI}
