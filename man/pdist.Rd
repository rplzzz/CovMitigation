% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run-mc.R
\name{pdist}
\alias{pdist}
\alias{genbeta}
\alias{gengamma}
\alias{genpois}
\alias{gennorm}
\title{Helper functions for creating parameter distributions}
\usage{
genbeta(alpha, beta)

gengamma(shape, rate)

genpois(lambda)

gennorm(mean, sd)
}
\arguments{
\item{alpha, beta}{Beta distribution parameters}

\item{shape, rate}{Gamma distribution parameters}

\item{lambda}{Poisson distribution parameter}

\item{mean, sd}{Normal distribution parameters}
}
\description{
These functions create parameter distributions for use in \code{\link{run_mc}}.
}
\section{Functions}{
\itemize{
\item \code{genbeta}: Create a beta distribution with specified parameters.

\item \code{gengamma}: Create a gamma distribution with specified shape and rate.

\item \code{genpois}: Create a Poisson distribution with specified mean.

\item \code{gennorm}: Create a normal distribution with specified mean and standard deviation.
}}

