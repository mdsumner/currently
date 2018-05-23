set.seed(7)
num_particles <- 150L
library(raster)
U <- raster("/rdsi/PRIVATE/raad/data_local/aad.gov.au/currents/fileu_00001.grd")
V <- raster("/rdsi/PRIVATE/raad/data_local/aad.gov.au/currents/filev_00001.grd")
U[is.na(U)] <- 0
V[is.na(V)] <- 0
n <- dim(U)[1L]
m <- dim(U)[2L]
particles <- cbind(x = runif(num_particles,min=0.501,max=(n + 0.499)),
                   y = runif(num_particles,min=0.501,max=(m + 0.499)))

t_step <- 1/(10*mean(c(abs(values(U)),abs(values(V)))))

xy <- list(x = particles[,1], y = particles[,2])
plot(flip(setExtent(sqrt(U*U + V* V), extent(0, ncol(U), 0, nrow(U))), "y"), col = viridis::viridis(32))
plot(xy, pch = ".")
vU <- values(U)
vV <- values(V)
L <- vector("list", 1e6/30)
for (i in 1:1e6) {
  xy <- rk:::rk(xy$x, xy$y,
                as.double(vU), as.double(vV),
                as.integer(m), as.integer(n), as.double(t_step))
  if (i %% 30 == 0) L[[i%/%30]] <- xy

}
