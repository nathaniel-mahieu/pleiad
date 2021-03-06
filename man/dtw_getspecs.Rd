% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getspecs.R
\name{dtw_getspecs}
\alias{dtw_getspecs}
\title{Get spectra after DTW alignment}
\usage{
dtw_getspecs(ms_list, align.mz, align.trace, align.rtrange, align.rtuserange,
  window.size)
}
\arguments{
\item{ms_list}{Raw data list containing list items k and s.}

\item{align.mz}{The mass in ms_list to align to.}

\item{align.trace}{Trace data.table containing columns rt and i to align to.}

\item{align.rtrange}{Vector of length 2 with begenning and end of retention time range to return. Range referse to the trace in align.trace.}
}
\value{
List containing specs, rt_map, and params.
}
\description{
\code{dtw_getspecs} aligns a mass trace to a supplied trace and returns a specified retention time range post alignment.
}
