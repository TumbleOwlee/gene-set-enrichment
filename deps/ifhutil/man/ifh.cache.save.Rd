% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util.R
\name{ifh.cache.save}
\alias{ifh.cache.save}
\title{Save the given variables in a cache file}
\usage{
ifh.cache.save(..., file)
}
\arguments{
\item{...}{arguments passed down to save(..)}

\item{file}{output file name}
}
\description{
This functions saves all given variables in a cache file. Currently it
is a simple wrapper around save().
}
