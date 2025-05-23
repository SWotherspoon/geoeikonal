% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geoeikonal.R
\name{geoeikonal}
\alias{geoeikonal}
\title{Solve the Eikonal Equation on a Lon/Lat Grid}
\usage{
geoeikonal(cost, extent)
}
\arguments{
\item{cost}{A matrix of travel costs - see details.}

\item{extent}{A vector of grid extents}
}
\value{
A matrix of minimum travel costs to a source point
}
\description{
Solve the Eikonal Equation on a Lon/Lat Grid.
}
\details{
This functions computes the minimum travel cost to each point in a
regular longitude/latitude grid from a set of source points by
solving the Eikonal equation with a fast marching method.

The \code{cost} argument is a matrix of travel costs through the cells
of the grid, and is oriented so that longitude increases down the
rows of the matrix, and latitude decreases across the columns.
Positive entries of the \code{cost} matrix give the cost of travel
through the corresponding cell of the grid, while \code{0} entries
represent source points, and negative, \code{NA} or \code{Inf} entries
represent barriers that cannot be traversed. The \code{extent} argument
is a vector of the form \code{c(lon_min, lon_max, lat_min, lat_max)}
that defines the geographic extent of the grid, assuming the
entries of the \code{cost} matrix represent the centres of the grid
cells.  form \code{c(lon_min, lon_max, lat_min, lat_max)} that defines
the geographic extent of the grid, assuming the entries of the
\code{cost} matrix represent the centres of the grid cells.

Distances are calculated assuming a spherical earth. If the
cost in each cell is one, then the results are distances in
kilometres.
}
\references{
Sethian, J.A. (1996), A fast marching level set method
for monotonically advancing fronts,
\emph{Proc. Natl. Acad. Sci.}  93 (4), 1591-1595.
}
