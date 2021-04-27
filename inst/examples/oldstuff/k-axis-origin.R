library(raadtools)
library(progress)
library(dplyr)

vlen <- function(x) sqrt(x[[1]]^2 + x[[2]]^2)
pvec <- function(x) plot(vlen(x), col = viridis::viridis(100))


## files
fds <- "/rdsi/PRIVATE/raad/data_local/aad.gov.au/currents/polar"
uregx <- "dt_south_polar_u_.*grd"
vregx <- "dt_south_polar_v_.*grd"
ufiles <- fs::dir_ls(fds, recursive = TRUE, type = "file", regexp = uregx)
vfiles <- fs::dir_ls(fds, recursive = TRUE, type = "file", regexp = vregx)
uvdates <- as.POSIXct(strptime(basename(ufiles), "dt_south_polar_u_%Y%m%d"))




## params
t_step <- 12 * 3600
reset <- 24 * 3600
startdate <- as.POSIXct("2016-02-20")
simulation_date <- startdate

## starter
ifile <- findInterval(startdate, uvdates)
read_current <- function(i) {
  disaggregate(brick(raster(ufiles[i]), raster(vfiles[ifile])), fact = 4)
}
curr <- read_current(ifile)

v <- tibble::as_tibble(readRDS("kaxis.rds"))
v$year <- 2016
v$date <- as.POSIXct(ISOdate(v$year, v$mon, v$day), tz = "UTC")
v[c("x", "y")] <- rgdal::project(as.matrix(v[c("lon", "lat")]), projection(curr))

## we want equal-distance approximate points along the entire track
v$distance <- c(0, trip::trackDistance(as.matrix(v[c("lon", "lat")])))
n_steps <- round(sum(v$distance)/(1.852 * 20))


avoy <- tibble(x = approxfun(cumsum(v$distance), v$x)(seq(0, sum(v$distance), length.out = n_steps)),
               y = approxfun(cumsum(v$distance), v$y)(seq(0, sum(v$distance), length.out = n_steps)),
               date = ISOdatetime(1970, 1, 1, 0, 0, 0, tz = "UTC") +
                 approxfun(cumsum(v$distance), v$date)(seq(0, max(v$distance), length.out = n_steps)))




nsamples <- 80
vpts <- avoy[c("x", "y", "date")] %>% mutate(track_i = row_number()) %>% slice(rep(track_i, nsamples)) %>% arrange(track_i)

vpts$x <- jitter(vpts$x, mean(range(vpts$x) * 5e-4))
vpts$y <- jitter(vpts$y, mean(range(vpts$y) * 5e-4))
plot(vpts$x, vpts$y, pch = ".")


vpts$sample <- rep(1:nrow(avoy), nsamples)
xy <- as.matrix(vpts[c("x", "y")])
plot(xy, asp = 1)

nn <- 141 * 2
l <- vector('list', nn)
i <- 0
ifile <- findInterval(startdate, uvdates)
method <- "simple"
elapsed <- 0

## need to offset each day's time

## keep proper simulation date

xy0 <- xy
pb <- progress::progress_bar$new(total = nn)
while(i < nn) {
  i <- i + 1
  ## first coefficient
  uv0 <- extract(curr, xy, method = method)
  ## second coefficient
  uv1 <- extract(curr, xy - uv0 * t_step / 2, method = method)
  ## third coefficient
  uv2 <- extract(curr, xy - uv1 * t_step / 2, method = method)
  uv3 <- extract(curr, xy - uv2 * t_step, method = method)
  pti <- xy + (uv0/6 + uv1/3 + uv2/3 + uv3/6) * t_step
  xy <- pti
  simulation_date <- simulation_date - t_step
  ## zap any locations that are not yet in time range
  too_early <- vpts$date < simulation_date
  xy[too_early, ] <- xy0[too_early, ]
  l[[i]] <- tibble(x = xy[,1], y = xy[,2], date = simulation_date, track_i = vpts$track_i, sample = vpts$sample)

  elapsed <- elapsed + t_step
  ifile <- findInterval(simulation_date, uvdates)
  curr <-  read_current(ifile)
  pb$tick()
}


x <- bind_rows(l,  .id = "iter")
#saveRDS(x, "particle_cache.rds")
library(ggplot2)
cl <- spTransform(rasterToContour(SOmap::Bathy, level = c(-4000, -2000, -1000, -5000)), projection(curr))
cltab <- fortify(spTransform(crop(cl, extent(range(v$x), range(v$y)) + 5e5), "+init=epsg:4326"))
x[c("lon", "lat")] <- rgdal::project(as.matrix(x[c("x", "y")]), projection(curr), inv = TRUE)
ggplot(x %>% sample_frac(0.25), aes(lon, lat, group = paste(track_i, sample), colour = date)) + geom_path() +
  geom_path(data = v, aes(lon, lat, group = NULL), colour = "black") + guides(colour = F) +
  geom_path(data = cltab, aes(long, lat, group = group), colour = "firebrick") +
  xlim(range(v$lon) + c(-1, 1) * 8) + ylim(range(v$lat) + c(-1, 1) * 5)







