% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transtwd97_function.R
\name{WGS84toTWD97}
\alias{WGS84toTWD97}
\title{Transforming your WGS84 data to TWD97 data}
\usage{
WGS84toTWD97(inputLat = NULL, inputLng = NULL)
}
\arguments{
\item{inputLat}{A numeric vector, WGS84 latitude of yours.}

\item{inputLng}{A numeric vector, WGS84 longitude of yours.}
}
\description{
This function allows you to transform your WGS84 location data to TWD97 data and create a new data frame.
WARNING: ONLY FOR LOCATION OF TAIWAN!!!!!!
}
\examples{
newTWD97data <- WGS84toTWD97(inputLat = 24.321, inputLng = 121.234)
}
\keyword{WGS84toTWD97}
