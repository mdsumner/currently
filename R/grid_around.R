
#' Grid around
#'
#' create a grid of points around a point (or points) in a matrix long,lat.
#'
#' @param x a point or points, 'cbind(longitude, latitude)'
#' @param radius width to buffer around the points (in degrees)
#' @param nxy size of grid 'c(nx, ny)'
#'
#' @return matrix of points longitude,latitude
#' @export
#'
#' @examples
#' grid_around(cbind(0, 0))
grid_around <- function(x, radius = 1, nxy = c(6, 6)) {
  xx <- range(x[,1])
  yy <- range(x[,2])
  raster::coordinates(raster::raster(raster::extent(xx[1] - radius, xx[2] + radius,
                                                    yy[1] - radius, yy[2] + radius),
                                     nrows = nxy[2], ncols = nxy[1]))
}
