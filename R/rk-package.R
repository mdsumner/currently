#' @keywords internal
"_PACKAGE"

#' @useDynLib 'rk', .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Particle tracking arena
#'
#' Example data from 'BioMassTracking' package.
#'
#' Comprised of a list with components
#' \itemize{
#' \item[lon] longitude vector
#' \item[lat] latitude vector
#' \item[U] matrix of U velocity component
#' \item[V] matrix of V velocity component
#' \item[S] matrix of S identifiers of regions
#' }
#' Orientation of the matrices is transpose to the `image` convention.
#' @docType data
#' @name arena
#' @examples
#' image(arena$lon, arena$lat, t(sqrt(arena$U ^2 + arena$V ^2)))
#' title("Kergulen Plateau, magnitude of surface currents")
NULL


#' Runge Kutta for particle tracking
#'
#' @param x x
#' @param ... dots
#'
#' @return 1
#' @export
#'
#' @examples
#' rk()
rk0 <- function(x, ...) {
  1
}

rk_run1 <- function(xy, ...) {

 # cbind(rk[[1]], rk[[2]])
  1
}
