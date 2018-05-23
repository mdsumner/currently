#library(raadtools)
#curr <- readAll(readcurr(xylim = extent(100, 150, -60, -45)))

ufiles <- list.files("/rdsi/PRIVATE/raad/data_local/aad.gov.au/currents", pattern = "fileu.*grd", full.names = TRUE)
vfiles <- list.files("/rdsi/PRIVATE/raad/data_local/aad.gov.au/currents", pattern = "filev.*grd", full.names = TRUE)

i <- 1
library(raster)
curr <- crop(brick(lapply(c(ufiles[1], vfiles[1]), raster)), e)

dur <- 24 * 3600

#curr <- crop(curr, e)
validcells <- which(!is.na(values(curr[[1]])))
cell <- sort(sample(validcells, 1000))
#cell <- cell[!is.na(extract(curr[[1]], cell))]
vec <- function(x) plot(sqrt(x[[1]]^2 + x[[2]]^2), col = viridis::viridis(100))

rk <- T
xy <- xyFromCell(curr, cell)
plot(xy, pch = ".", asp = 1)

t_step <- 1/(1*mean(c(abs(values(curr[[1]])),
                      abs(values(curr[[2]]))), na.rm = TRUE))

t_step <- 24 * 3600 * 10
i <- 0
lastxy <- xy
totaldur <- 0
while(TRUE) {

  ## first coefficient
  uv0 <- extract(curr, cell)
  if (rk) {
    ## second coefficient
    uv1 <- extract(curr, xy + uv0 * t_step / 2, method = "bilinear")
    ## third coefficient
    uv2 <- extract(curr, xy + (uv1) * t_step / 2, method = "bilinear")
    uv3 <- extract(curr, xy + uv2 * t_step, method = "bilinear")
    pti <- xy + uv0/6 + uv1/3 + uv2/3 + uv3/6
  } else {
    pti <- xy + uv0
  }
  if ((totaldur %% dur) < 1200) {
    i <- i + 1
    curr <- crop(brick(lapply(c(ufiles[i], vfiles[i]), raster)), e)
    print(i)

    segments(lastxy[,1], lastxy[,2], pti[,1], pti[,2])
    bad <- is.na(xy[,1])
    if (any(bad)) xy[bad, ] <- xyFromCell(curr, sample(validcells, sum(bad)))
    lastxy <- xy
  }
  totaldur <- totaldur + t_step/2400

  xy <- pti
}
