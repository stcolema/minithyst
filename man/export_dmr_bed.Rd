% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmr_utils.R
\name{export_dmr_bed}
\alias{export_dmr_bed}
\title{Export DMRs to BED format}
\usage{
export_dmr_bed(
  dmr_results,
  output_file,
  track_name = "DMRs",
  separate_by = NULL
)
}
\arguments{
\item{dmr_results}{Filtered or collapsed DMR results}

\item{output_file}{Path to output BED file}

\item{track_name}{Name for the track}

\item{separate_by}{Character, column to separate tracks by (e.g., "direction")}
}
\description{
Export DMRs to BED format
}
