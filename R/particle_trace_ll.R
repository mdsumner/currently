xyuv2l <- function(xy, x) {
  xscale <- 1/cos(xy[,2] * pi/180)

  cbind(x[,1] * xscale, x[,2]) / (1852 * 60)

}

#' @name particle_trace
#' @export
particle_trace_ll <- function(xy,
                           time_step = 24 * 3600,
                           start_date = NULL,
                           end_date = NULL,
                           method = "bilinear", plot = FALSE, silent = FALSE, rk = FALSE) {

  readcurr0 <- memoise::memoize(raadtools:::readcurr)

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
  files <- raadfiles::altimetry_daily_files()
  for (jj in seq_along(dates)) {
    model_time <- dates[jj]
    curr <- readcurr0(model_time, inputfiles = files)
    # if (jj == 1L) {
    #   xy <-   reproj::reproj_xy(xy, raster::projection(curr), source = "+proj=longlat +datum=WGS84")
    #   if (plot) plot(xy, asp = 1)
    # }
    #
    ##i <- i + 1
    ## first coefficient

    ## convert uv metres to longlat

    uv0 <- raster::extract(curr, xy, method = method)
    if (rk) {
      ## second coefficient
      uv1 <- extract(curr, xy + (uv0 * time_step)/2, method = method)
      ## third coefficient
      uv2 <- extract(curr, xy + (uv1 * time_step)/2, method = method)
      uv3 <- extract(curr, xy + (uv2 * time_step), method = method)
      uv1 <- xyuv2l(xy, uv1)
      uv2 <- xyuv2l(xy, uv2)
      uv3 <- xyuv2l(xy, uv2)

      pti <- xy + (uv0/6 + uv1/3 + uv2/3 + uv3/6) * time_step
    } else {
      pti <- xy + xyuv2l(xy, uv0 )* time_step
    }
    xy <- pti
    l[[jj]] <- xy

    if (!silent) pb$tick()
    if (plot) points(xy, pch = ".")
  }

  pts <- setNames(as.data.frame(do.call(rbind, l)), c("lon", "lat"))
  pts$group <- rep(1:nrow(l[[1]]), length.out = nrow(pts))
  pts$date <- rep(dates, lengths(l)/ncol(l[[1]]))
  tibble::as_tibble(pts)
}
