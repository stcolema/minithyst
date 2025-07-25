% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{diagnose_index_metadata}
\alias{diagnose_index_metadata}
\title{Diagnose index and metadata compatibility}
\usage{
diagnose_index_metadata(chr_index, metadata, h5_paths)
}
\arguments{
\item{chr_index}{Optimized index object}

\item{metadata}{data.table with cell metadata}

\item{h5_paths}{data.table with h5 file paths}
}
\value{
Invisible list with diagnostic information
}
\description{
Checks for common issues with cell ID matching between index, metadata,
and h5_paths. Useful for debugging when calc_smoothed_windows fails.
}
\examples{
\dontrun{
# Run diagnostic if calc_smoothed_windows fails
diagnose_index_metadata(index, metadata, h5_paths)
}
}
