#' Particle tracing
#'
#' Create traces of particles from velocity vector fields. This assumes availability of the raadtools package.
#'
#' Relies on polar versions of Copernicus altimetry surface velocity fields.
#'
#' This function is the high-level vectorized "use the raster package" version of some older C-based particle
#' tracing code. The technique uses raster to extract the x and y displacements in surface current vectors, and
#' apply those to input points at the given time step. Runge-Kutta can be used to generate 4 coefficients for the
#' displacements.
#'
#' We can use 'start_date' and 'end_date' to specify the direction, the start is the time we begin the
#' simulation at, and the end when it finishes. If the 'end_date' is prior to the 'start_date' the time step
#' is negated internally.
#'
#' @param xy input point/s in longitude,latitude
#' @param time_step duration of model time step (positive, in seconds)
#' @param start_date time to start the model (can be of type Date but beware doing your own math mixing Date/POSIXct)
#' @param end_date time to end the model (may be prior to or preceding start_date)
#' @param method in the raster lookup, may be 'bilinear' or 'simple'
#' @param plot plot the points or not (makes it slower)
#' @param silent suppress any message
#' @param rk use multiple coefficient Runge-Kutta (default is `FALSE`)
#'
#' @importFrom memoise memoize
#' @importFrom progress progress_bar
#' @importFrom raster extract projection
#' @importFrom rgdal project
#' @importFrom stats setNames
#' @importFrom tibble as_tibble
#' @importFrom graphics plot points
#' @export
particle_trace <- function(xy,
                           time_step = 24 * 3600,
                           start_date = NULL,
                           end_date = NULL,
                           method = "bilinear", plot = FALSE, silent = FALSE, rk = FALSE) {

  readcurr_POLAR <- memoise::memoize(raadtools:::readcurr_polar)

  if (time_step < 0) stop("'time_step' must be positive (use 'start_date' and 'end_date' to specify direction)")

  if (is.null(start_date)) stop("must specify start_date e.g. as.Date(\"2021-02-15\")")
  start_date <- as.POSIXct(start_date, tz = "UTC")
  if (is.null(end_date)) stop("must specify end_date e.g. as.Date(\"2021-02-15\")")
  end_date <- as.POSIXct(end_date, tz = "UTC")
  if (end_date == start_date) stop("'end_date' must not be equal to 'start_date'")
  if (end_date < start_date)  time_step <- -time_step


  ## these are in forward or reverse time depending on the relation of start to end
  dates <- seq(start_date, end_date, by = time_step)

  N <- length(dates)

  if (!silent) {
    message(sprintf("proceeding with %s time step from %s to %s with %i increments",
                    c("positive", "negative")[(sign(time_step) < 0) + 1],
                    as.character(start_date), as.character(end_date), N))
  }

  l <- vector("list", N)
  ## ticker
  if (!silent) pb <- progress::progress_bar$new(total = N)
  files <- raadfiles::altimetry_currents_polar_files()
  for (jj in seq_along(dates)) {
    model_time <- dates[jj]
    curr <- readcurr_POLAR(model_time, inputfiles = files)
    if (jj == 1L) {
      xy <-   rgdal::project(xy, raster::projection(curr))
      if (plot) plot(xy, asp = 1)
    }

    ##i <- i + 1
    ## first coefficient
    uv0 <- extract(curr, xy, method = method)
    if (rk) {
      ## second coefficient
      uv1 <- extract(curr, xy + (uv0 * time_step)/2, method = method)
      ## third coefficient
      uv2 <- extract(curr, xy + (uv1 * time_step)/2, method = method)
      uv3 <- extract(curr, xy + (uv2 * time_step), method = method)
      pti <- xy + (uv0/6 + uv1/3 + uv2/3 + uv3/6) * time_step
    } else {
      pti <- xy + uv0 * time_step
    }
    xy <- pti
    l[[jj]] <- xy

    if (!silent) pb$tick()
    if (plot) points(xy, pch = ".")
  }

  pts <- setNames(tibble::as_tibble(do.call(rbind, l)), c("x", "y"))
  pts$group <- rep(1:nrow(l[[1]]), length.out = nrow(pts))
  pts[c("lon", "lat")] <- rgdal::project(cbind(pts$x, pts$y), projection(curr), inv = T)
  pts$date <- rep(dates, lengths(l)/ncol(l[[1]]))
  pts
}
