% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SAS_shapes.R
\name{shape_step}
\alias{shape_step}
\title{step SAS shape}
\usage{
shape_step(P_s, parameters, ...)
}
\arguments{
\item{P_s}{the cumulative storage age distribution (S_T/S)}

\item{parameters}{list with parameter 'u', which indicates the top interval of
where Omega is linear. 'u' must be between 0 and 1.}

\item{...}{unused arguments}
}
\value{
Returns the cumulative SAS function, Omega.
}
\description{
A function to compute a step fSAS function (a "left" fSAS function).
}
\details{
When P_s < u, Omega = P_u/u. When P_s >= u, Omega = 1.
}
