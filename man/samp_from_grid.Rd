% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samp_from_grid.R
\name{samp_from_grid}
\alias{samp_from_grid}
\title{Sample from grid approximation of posterior distribuitions}
\usage{
samp_from_grid(joint_post, params, nsamps)
}
\arguments{
\item{joint_post}{Data frame giving the joint posterior}

\item{params}{Vector of names of parameters to sample}

\item{nsamps}{Number of samples to take}
}
\value{
A data.frame that is nsamps by length(params) of samples
}
\description{
Sample from the joint posterior of parameters
}
