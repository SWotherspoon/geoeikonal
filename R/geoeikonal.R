##' Solve the Eikonal Equation on a Lon/Lat Grid.
##'
##' This functions computes the minimum travel cost to each point in a
##' regular longitude/latitude grid from a set of source points by
##' solving the Eikonal equation with a fast marching method.
##'
##' The `cost` argument is a matrix of travel costs through the cells
##' of the grid, and is oriented so that longitude increases down the
##' rows of the matrix, and latitude decreases across the columns.
##' Positive entries of the `cost` matrix give the cost of travel
##' through the corresponding cell of the grid, while `0` entries
##' represent source points, and negative, `NA` or `Inf` entries
##' represent barriers that cannot be traversed. The `extent` argument
##' is a vector of the form `c(lon_min, lon_max, lat_min, lat_max)`
##' that defines the geographic extent of the grid, assuming the
##' entries of the `cost` matrix represent the centres of the grid
##' cells.  form `c(lon_min, lon_max, lat_min, lat_max)` that defines
##' the geographic extent of the grid, assuming the entries of the
##' `cost` matrix represent the centres of the grid cells.
##' 
##' Distances are calculated assuming a spherical earth. If the 
##' cost in each cell is one, then the results are distances in 
##' kilometres.
##'
##' @title Solve the Eikonal Equation on a Lon/Lat Grid
##' @param cost A matrix of travel costs - see details.
##' @param extent A vector of grid extents
##' @return A matrix of minimum travel costs to a source point
##' @references Sethian, J.A. (1996), A fast marching level set method
##'   for monotonically advancing fronts,
##'   \emph{Proc. Natl. Acad. Sci.}  93 (4), 1591-1595.
##' @export
geoeikonal <- function(cost, extent) {
  ## Validate inputs
  if(!is.matrix(cost)) stop("cost must be a matrix")
  if(!is.vector(extent) || length(extent) != 4 ||
       extent[1] >= extent[2] || extent[3] >= extent[4] ||
       extent[3] < -90 || extent[4] > 90)
    stop("extent must be a vector of the form c(lon_min, lon_max, lat_min, lat_max)")
  if(!is.double(cost)) storage.mode(cost) <- "double"
  
  .Call(C_geoeikonal, cost, extent, matrix(0,nrow(cost),ncol(cost)))
}



##' Solve the Eikonal Equation on a regular Grid.
##'
##' This functions computes the minimum travel cost to each point in a
##' regular grid from a set of source points by solving the Eikonal
##' equation with a fast marching method.
##'
##' The `cost` argument is a matrix of travel costs through the cells
##' of the grid, and is oriented so that x increases down the rows of
##' the grid, and y increases across the columns.  Positive entries of
##' the `cost` matrix give the cost of travel through the
##' corresponding cell of the grid, while `0` entries represent a
##' source point, and negative, `NA` or `Inf` entries represent
##' barriers that cannot be traversed.  The `hx` and `hy` arguments
##' give the grid spacings in the x and y directions.
##'
##' @title Solve the Eikonal Equation on a Regular Grid
##' @param cost A matrix of travel costs - see details.
##' @param hx The grid spacing in the x direction
##' @param hy The grid spacing in the y direction
##' @return A matrix of minimum travel costs to a source point
##' @references Sethian, J.A. (1996), A fast marching level set method
##'   for monotonically advancing fronts,
##'   \emph{Proc. Natl. Acad. Sci.}  93 (4), 1591-1595.
##' @export
grdeikonal <- function(cost, hx=1, hy=1) {
  ## Validate inputs
  if(!is.matrix(cost)) stop("cost must be a matrix")
  if(!is.double(cost)) storage.mode(cost) <- "double"
  
  .Call(C_grdeikonal, cost, hx, hy, matrix(0.0,nrow(cost),ncol(cost)))
}


## Solve the Eikonal equation for a lon/lat raster.
##'
##' This functions computes the minimum travel cost to each point in a
##' longitude/latitude raster from a set of source points by solving
##' the Eikonal equation with a fast marching method.
##'
##' The `cost` argument must be a `SpatRaster` with a
##' longitude/latitude coordinate reference system.  Positive entries
##' of the `cost` raster give the cost of travel through the
##' corresponding cell of the raster, while `0` entries represent a
##' source point, and negative, `NA` or `Inf` entries represent
##' barriers that cannot be traversed.
##'
##' Distances are calculated assuming a spherical earth. If the 
##' cost in each cell is one, then the results are distances in 
##' kilometres.
##'
##' @title Solve the Eikonal Equation on a Lon/Lat Raster
##' @param cost A `SpatRaster` with a longitude/latitude coordinate
##'   reference system.
##' @return A matrix of minimum travel costs to a source point
##' @export
rsteikonal <- function(cost) {
  ## Check CRS is long/lat and get extent
  if(!terra::is.lonlat(cost)) stop("cost matrix must have a lon/lat CRS")
  extent <- as.vector(terra::ext(cost))
  ## Convert the raster to a matrix
  C <- matrix(terra::values(cost),terra::ncol(cost),terra::nrow(cost))
  r <- .Call(C_geoeikonal, C, extent, matrix(0,nrow(C),ncol(C)))
  ## Return as raster
  terra::rast(nrows=terra::nrow(cost), ncols=terra::ncol(cost),
              extent=terra::ext(cost), crs=terra::crs(cost), 
              vals=as.vector(r))
}


##' Compute the (upwind) gradient of a solution to the Eikonal equation 
##' computed with [geoeikonal()].
##'
##' @title Upwind Gradient of Distance on a Lon/Lat Grid
##' @param dist A matrix of distances computed with [geoeikonal()].
##' @param extent A vector of grid extents
##' @return A two element list containing the components of the gradient
##' @export
geogradient <- function(dist, extent) {
  ## Validate inputs
  if(!is.matrix(dist)) stop("dist must be a matrix")
  if(!is.vector(extent) || length(extent) != 4 ||
       extent[1] >= extent[2] || extent[3] >= extent[4] ||
       extent[3] < -90 || extent[4] > 90)
    stop("extent must be a vector of the form c(lon_min, lon_max, lat_min, lat_max)")
  if(!is.double(dist)) storage.mode(dist) <- "double"
  
  .Call(C_geogradient, dist, extent, matrix(0,nrow(dist),ncol(dist)), matrix(0,nrow(dist),ncol(dist)))
}


##' Compute the (upwind) gradient of a solution to the Eikonal equation 
##' computed with [grdeikonal()].
##'
##' @title Upwind Gradient of Distance on a Lon/Lat Grid
##' @param dist A matrix of distances computed with [grdeikonal()].
##' @param hx The grid spacing in the x direction
##' @param hy The grid spacing in the y direction
##' @return A two element list containing the components of the gradient
##' @export
grdgradient <- function(dist, hx=1, hy=1) {
  ## Validate inputs
  if(!is.matrix(dist)) stop("dist must be a matrix")
  if(!is.double(dist)) storage.mode(dist) <- "double"
  
  .Call(C_grdgradient, dist, hx, hy, matrix(0,nrow(dist),ncol(dist)), matrix(0,nrow(dist),ncol(dist)))
}

##' Compute the (upwind) gradient of a solution to the Eikonal equation 
##' computed with [rsteikonal()].
##'
##' @title Upwind Gradient of Distance on a Lon/Lat Raster
##' @param dist A matrix of distances computed with [rsteikonal()].
##' @return A two element list containing the components of the gradient
##' @export
rstgradient <- function(dist) {
  ## Check CRS is long/lat and get extent
  if(!terra::is.lonlat(dist)) stop("dist matrix must have a lon/lat CRS")
  extent <- as.vector(terra::ext(dist))
  ## Convert the raster to a matrix
  D <- matrix(terra::values(dist),terra::ncol(dist),terra::nrow(dist))
  
  r <-.Call(C_geogradient, D, extent, matrix(0,nrow(D),ncol(D)), matrix(0,nrow(D),ncol(D)))
  terra::rast(nrows=terra::nrow(dist), ncols=terra::ncol(dist), nlyrs=2,
              extent=terra::ext(dist), crs=terra::crs(dist), 
              vals=as.vector(c(r[[1]],r[[2]])),names=c("lon","lat"))
}


