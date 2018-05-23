#tools::package_native_routine_registration_skeleton("../rk", "src/init.c",character_only = FALSE)


set.seed(7)
num_particles <- 500L
m <- dim(arena$U)[1L]
n <- dim(arena$V)[2L]
particles <- cbind(x = runif(num_particles,min=0.501,max=(n + 0.499)),
                   y = runif(num_particles,min=0.501,max=(m + 0.499)))
t_step <- 1/(10*mean(c(abs(arena$U),abs(arena$V))))
plot(particles)
for (i in 1:1e3) {
  xy <- rk:::rk(as.double(particles[,1L]), as.double(particles[,2L]),
           as.double(arena$U), as.double(arena$V),
           as.integer(m), as.integer(n), as.double(t_step))
  points(xy, pch = ".")
  particles <- do.call(cbind, xy)
  bad <- particles[,1] < 0.5 | particles[,1] >  n | particles[,2] < 0.5 | particles[,2] > m
  if (any(bad)) {
    nnew <- sum(bad)

    particles[bad, ] <- cbind(x = runif(nnew,min=0.501,max=(n + 0.499)),
                              y = runif(nnew,min=0.501,max=(m + 0.499)))
  }
  #if (i %% 1000 == 0) plot(particles)
}


## if U and V are rasters in longlat (southern hemisphere), then
arena <- prepare.arena(as.data.frame(U, xy = T) %>% arrange(y, x),
                       as.data.frame(V, xy = T) %>% arrange(y, x))

