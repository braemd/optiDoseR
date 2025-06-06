% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nm_update_model.R
\name{nmUpdateData}
\alias{nmUpdateData}
\title{Update $DATA in NONMEM model}
\usage{
nmUpdateData(optiProject, name)
}
\arguments{
\item{optiProject}{A optiProject-object with a model added.}

\item{new_data_name}{Name of the new OptiProject-data generated through
\code{DataGen}.}
}
\description{
Update the $DATA line in a NONMEM model to new optiProject data
}
\details{
Should be called after \code{nmUpdateInputs}
}
\author{
Dominic Bräm
}
