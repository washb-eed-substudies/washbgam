% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_RE_gam.R
\name{fit_RE_gam}
\alias{fit_RE_gam}
\title{Title}
\usage{
fit_RE_gam(
  d,
  Y,
  X,
  W = NULL,
  forcedW = NULL,
  V = NULL,
  id = "clusterid",
  family = "gaussian",
  pval = 0.2,
  print = TRUE
)
}
\arguments{
\item{forcedW}{Set to NA to skip forcing any variables, don't include to skip prescreening of "tr" and variables starting with age_, ageday_, or agedays_. Or manually set with a vector of variable names.}

\item{print}{}
}
\value{

}
\description{
Title
}
