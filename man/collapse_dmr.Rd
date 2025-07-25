% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmr_analysis.R
\name{collapse_dmr}
\alias{collapse_dmr}
\title{Collapse DMRs}
\usage{
collapse_dmr(dmr_filtered, max_gap = 0, min_length = 1000)
}
\arguments{
\item{dmr_filtered}{Filtered results from filter_dmr}

\item{max_gap}{Maximum gap between DMRs to merge}

\item{min_length}{Minimum length of collapsed DMR}
}
\value{
data.table with collapsed DMRs
}
\description{
Collapse adjacent DMRs into a single entity
}
