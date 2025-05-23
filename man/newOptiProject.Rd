% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adder.R
\name{newOptiProject}
\alias{newOptiProject}
\title{Create new optiProject}
\usage{
newOptiProject(software)
}
\arguments{
\item{software}{Software for which the optiProject should be initialized ("Monolix", "NONMEM", and "nlmixr2" available)}
}
\value{
A optiProject-object, to which additional elements of the optiDose approach
can be added.
}
\description{
Initializes a new optiProject-object for a software to which a model, parameters,
dose levels, and constraints/targets can be added.
}
\examples{
\dontrun{
warf_proj <- newOptiProject("Monolix")
}
}
\author{
Dominic Bräm
}
