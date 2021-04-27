library(raadtools)
ii <- readRDS("k-axis-particles.rds")
cfiles <- raadfiles::altimetry_currents_polar_files()
curr <- raster::raster(cfiles$ufullname[1])

ii[c("lon", "lat")] <- reproj::reproj(ii[c("x", "y")], source = raster::projection(curr), target = 4326)[,1:2]
bad <- is.na(ii$lon) | is.na(ii$lat)
ii$ghrsst <- NA
ex <- raster::extent(range(ii$lon[!bad]), range(ii$lat[!bad])) + 2
library(raadtools)
ii$ghrsst[!bad] <- extract(readghrsst, ii[!bad, c("lon", "lat", "date")])
saveRDS(ii, "k-axis-ghrsst.rds")

dd <- readRDS("k-axis-ghrsst.rds")

library(raadtools)
bad <- is.na(dd$lon) | is.na(dd$lat)
dd$ice <- rep(NA_real_, nrow(dd))
dd$ice[!bad] <- extract(readice, dd[!bad, c("lon", "lat", "date")])
dd$dist_ice_edge <- rep(NA_real_, nrow(dd))
dd$dist_ice_edge[!bad] <- extract(distance_to_ice_edge, dd[!bad, c("lon", "lat", "date")])


dd$DATE <- as.Date(dd$date)
ud <- unique(dd$DATE)
gg <- raster(extent(range(dd$lon, na.rm = T), range(dd$lat, na.rm = TRUE)), res = 0.5)
dd$chl <- rep(NA_real_, nrow(dd))
for (i in seq_along(ud)) {
  dts <- seq(ud[i] - 30, by = "1 day", length = 31)
  chl <- readchla(dts, grid = gg)
  isub <- dd$DATE == ud[i]
  xy <- dd[isub, c("lon", "lat")]
  bad <- is.na(dd[["lon"]]) | is.na(dd[["lat"]])
  if (!all(bad)) {
    dd$chl[isub][!bad] <- extract(focal(chl, matrix(1, 3, 3), fun = median, na.rm = T),
                          xy[!bad, ])
  }
  print(i)
}


dd$ghrsst <- dd$ghrsst - 273.15
saveRDS(dd, file = "k-axis-particles-sst-ice-chl.rds")
