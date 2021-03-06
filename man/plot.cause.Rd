% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.cause}
\alias{plot.cause}
\title{Plot CAUSE}
\usage{
\method{plot}{cause}(res, intern = FALSE, type = c("posteriors", "data"), pval_thresh = 5e-08)
}
\arguments{
\item{res}{object of class cause}

\item{intern}{If TRUE, function returns a list of grobs. Otherwise it plots.}

\item{type}{Either "posteriors" or "data". See details.}

\item{pval_thresh}{p-value threshold used if type=="data".}
}
\description{
Plot posteriors for sharing and causal models, display summary tables
}
\details{
If type == "posteriors", the function will plot the posterior distributions of the parameters
and display summary tables giving medians and credible intervals. If type == "data" the function
will plot the data thresholded on the trait 1 p-value if using pval_thresh.
}
