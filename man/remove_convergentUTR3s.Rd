% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/14.remove_convergentUTR3s.R
\name{remove_convergentUTR3s}
\alias{remove_convergentUTR3s}
\title{remove the converging candidates 3' UTRs LIKE UTR3___UTR3}
\usage{
remove_convergentUTR3s(x)
}
\arguments{
\item{x}{the distal 3UTR coverage}
}
\value{
the 3UTR coverage after removing the next 3UTR
}
\description{
some of the results is from connected two 3' UTRs. We want to remove them.
}
\details{
The algorithm need to be improved.
}
\author{
Jianhong Ou
}
\keyword{internal}
