% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Kim et al.'s Modified T-test.R
\name{modified.t.test}
\alias{modified.t.test}
\title{Kim et al.'s Modified T-test}
\usage{
modified.t.test(x, y, alternative = "two.sided")
}
\arguments{
\item{x}{a (non-empty) numeric vector of data values with some missing value(NA).}

\item{y}{a (non-empty) numeric vector of data values with some missing value(NA).}

\item{alternative}{a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".}
}
\value{
Results of test including p-value will be printed.
}
\description{
Performs Kim et al.'s Modified T-test on a partially matched sample in a form of R vector.
}
\examples{
# Generating Toy Examples
set.seed(123)
x <- rnorm(20)
x[sample(1:20, 3)] <- NA # Deliverately generating some missing values
y <- (rnorm(20) + 1)/3
y[sample(which(!is.na(x)), 4)] <- NA
modified.t.test(x, y)
modified.t.test(x, y, alternative = "greater")
modified.t.test(x, y, alternative = "less")

}
