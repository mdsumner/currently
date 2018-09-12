
library(raadtools)
##curr <- readAll(readcurr(xylim = extent(100, 150, -60, -45)))

fds <- "/rdsi/PRIVATE/raad/data_local/aad.gov.au/currents"
ufiles <- list.files(fds, pattern = "fileu.*grd", full.names = TRUE)
vfiles <- list.files(fds, pattern = "filev.*grd", full.names = TRUE)

ifile <- 1
curr <- brick(raster(ufiles[ifile]), raster(vfiles[ifile]))
vlen <- function(x) sqrt(x[[1]]^2 + x[[2]]^2)

valid <- which(!is.na(values(curr[[1]])))
sampler <- function(n = 8e4) {
  sample(valid, n)
}
pvec <- function(x) plot(vlen(x), col = viridis::viridis(100))


# t_step <- 1/(10*mean(c(abs(values(curr[[1]])),
#                         abs(values(curr[[2]]))), na.rm = TRUE))
t_step <- 24 * 3600
reset <- 24 * 3600
startdate <- as.Date("2016-01-20")

#t_step <- 30 * 60
i <- 0
xy <- xyFromCell(curr, sampler(n = 500))

#pt <- cbind(2984923, -1358795)
#xy <- do.call(rbind, replicate(20, jitter(pt, 0.05), simplify = F))
#points(xy, col = "white")
l <- vector('list', length(ufiles))
i <- 0
method <- "simple"
elapsed <- 0
while(i < length(ufiles)) {
  ##if (i %% 500 == 0) xy <- xyFromCell(curr, sampler(n = 500))
  i <- i + 1
  ## first coefficient
  uv0 <- extract(curr, xy, method = method)
  ## second coefficient
  uv1 <- extract(curr, xy + uv0 * t_step / 2, method = method)
  ## third coefficient
  uv2 <- extract(curr, xy + uv1 * t_step / 2, method = method)
  uv3 <- extract(curr, xy + uv2 * t_step, method = method)
  pti <- xy + (uv0/6 + uv1/3 + uv2/3 + uv3/6) * t_step
  #segments(xy[,1], xy[,2], pti[,1], pti[,2])
  xy <- pti
  l[[i]] <- xy
  elapsed <- elapsed + t_step
  if (elapsed >= reset) {
    ifile <- ifile + 1
    curr <- brick(raster(ufiles[i]), raster(vfiles[i]))
    print(i)
    elapsed <- 0
  }
}


l1 <- purrr::map_df(l, ~setNames(tibble::as_tibble(.x), c("X", "Y")) %>% mutate(id = row_number()), .id = "iter")

saveRDS(l1, sprintf("particles/particles_%04i.rds", j))
x <- l1 %>% mutate(iter = as.integer(iter)) %>% arrange(id, iter) %>% transmute(x_ = X, y_ = Y, order_ = row_number(), branch_ = id, object_ = id) %>%
  dplyr::filter(!is.na(x_), !is.na(y_)) %>% spbabel::sp()



library(DBI)
dbdir <- tempdir()
library(MonetDBLite)
con <- dbConnect(MonetDBLite::MonetDBLite(), dbdir)

## this works!
x <- tbl(con, "particle") %>% dplyr::filter(id < 100, !is.na(`X`))  %>% transmute(x_  = `X`, y_ = `Y`, order_ = row_number(), branch_ = 1, object_ = id) %>% collect() %>% spbabel::sp(crs = projection(curr))
dbWriteTable(con, "particle", l1)

