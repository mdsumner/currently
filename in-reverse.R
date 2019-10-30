library(raadtools)

files <- raadfiles::altimetry_currents_polar_files()
ifile <- 1

curr <- brick(raster(files$ufullname[ifile]), raster(files$vfullname[ifile]))
vlen <- function(x) sqrt(x[[1]]^2 + x[[2]]^2)

valid <- which(!is.na(values(curr[[1]])))
sampler <- function(n = 8e4) {
  sample(valid, n)
}
pvec <- function(x) plot(vlen(x), col = viridis::viridis(100))


# t_step <- 1/(10*mean(c(abs(values(curr[[1]])),
#                         abs(values(curr[[2]]))), na.rm = TRUE))
t_step <- 3600
reset <- 24 * 3600
startdate <- as.Date("2016-01-20")
ifile <- findInterval(startdate, as.Date(files$date))
#t_step <- 30 * 60
i <- 0
xy <- xyFromCell(curr, sampler(n = 100))

#pt <- cbind(2984923, -1358795)
#xy <- do.call(rbind, replicate(20, jitter(pt, 0.05), simplify = F))
plot(xy, col = "white")
nn <- 5000
l <- vector('list', nn)
i <- 0
method <- "simple"
elapsed <- 0
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
  segments(xy[,1], xy[,2], pti[,1], pti[,2])
  xy <- pti
  l[[i]] <- xy
  elapsed <- elapsed + t_step
  if (elapsed >= reset) {
    ifile <- ifile- 1
    brick(raster(files$ufullname[ifile]), raster(files$vfullname[ifile]))
    print(i)
    elapsed <- 0
  }
}


library(ggplot2)
library(dplyr)
purrr::map_df(l, ~setNames(tibble::as_tibble(.x), c("X", "Y")) %>% mutate(id = row_number()), .id = "iter") %>% ggplot() +
  aes(X, Y, group = id, colour = id) + geom_path()


