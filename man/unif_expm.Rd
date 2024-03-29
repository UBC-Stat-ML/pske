% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unif_expm.R
\name{unif_expm}
\alias{unif_expm}
\title{Uniformization}
\usage{
unif_expm(Q, K, t_pow = 1, eps, sparse)
}
\arguments{
\item{Q}{rate matrix.}

\item{K}{precision level. If provided, eps is ignored.}

\item{t_pow}{the t in exp(tQ)}

\item{eps}{desired error tolerance}
}
\value{
exp(tQ) up to tolerance eps in op-1 norm.
}
\description{
Uniformization
}
