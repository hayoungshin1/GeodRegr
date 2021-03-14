pkgname <- "GeodRegr"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('GeodRegr')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("are")
### * are

flush(stderr()); flush(stdout())

### Name: are
### Title: Approximate ARE of an M-type estimator to the least-squares
###   estimator
### Aliases: are

### ** Examples

are('l1', 10)




cleanEx()
nameEx("are_nr")
### * are_nr

flush(stderr()); flush(stdout())

### Name: are_nr
### Title: Newton-Raphson method for the 'are' function
### Aliases: are_nr

### ** Examples

dimension <- 4
x <- 1:10000 / 1000
# use a graph of the are function to pick a good starting point
plot(x, are('huber', dimension, x) - 0.95)
are_nr('huber', dimension, 2)




cleanEx()
nameEx("calvaria")
### * calvaria

flush(stderr()); flush(stdout())

### Name: calvaria
### Title: Data on calvaria growth in rat skulls
### Aliases: calvaria
### Keywords: datasets

### ** Examples

# we will test the robustness of each estimator by comparing their
# performance on the original (corrupted) data set to that of the L_2
# estimator on the uncorrupted data set (with the 4 problematic data points
# removed).

data(calvaria)

manifold <- 'kendall'

contam_x_data <- calvaria$x
contam_mean_x <- mean(contam_x_data)
contam_x_data <- contam_x_data - contam_mean_x # center x data
uncontam_x_data <- calvaria$x[ -c(23, 101, 104, 160)]
uncontam_mean_x <- mean(uncontam_x_data)
uncontam_x_data <- uncontam_x_data - uncontam_mean_x # center x data

contam_y_data <- calvaria$y
uncontam_y_data <- calvaria$y[, -c(23, 101, 104, 160)] # remove corrupted
    # columns

landmarks <- dim(contam_y_data)[1]
dimension <- 2 * landmarks - 4

# we ignore Huber's estimator as the L_1 estimator already has an
# (approximate) efficiency above 95% in 12 dimensions; see documentation for
# the are and are_nr functions

tol <- 1e-5
uncontam_l2 <- geo_reg(manifold, uncontam_x_data, uncontam_y_data,
    'l2', p_tol = tol, V_tol = tol)
contam_l2 <- geo_reg(manifold, contam_x_data, contam_y_data,
    'l2', p_tol = tol, V_tol = tol)
contam_l1 <- geo_reg(manifold, contam_x_data, contam_y_data,
    'l1', p_tol = tol, V_tol = tol)
contam_tukey <- geo_reg(manifold, contam_x_data, contam_y_data,
    'tukey', are_nr('tukey', dimension, 10, 0.99), p_tol = tol, V_tol = tol)

geodesics <- vector('list')
geodesics[[1]] <- uncontam_l2
geodesics[[2]] <- contam_l2
geodesics[[3]] <- contam_l1
geodesics[[4]] <- contam_tukey

loss(manifold, geodesics[[1]]$p, geodesics[[1]]$V, uncontam_x_data,
    uncontam_y_data, 'l2')
loss(manifold, geodesics[[2]]$p, geodesics[[2]]$V, contam_x_data,
    contam_y_data, 'l2')
loss(manifold, geodesics[[3]]$p, geodesics[[3]]$V, contam_x_data,
    contam_y_data, 'l1')
loss(manifold, geodesics[[4]]$p, geodesics[[4]]$V, contam_x_data,
    contam_y_data, 'tukey', are_nr('tukey', dimension, 10, 0.99))

# visualization of each geodesic

oldpar <- par(mfrow = c(1, 4))

days <- c(7, 14, 21, 30, 40, 60, 90, 150)
pal <- colorRampPalette(c("blue", "red"))(length(days))

# each predicted geodesic will be represented as a sequence of the predicted
# shapes at each of the above ages, the blue contour will show the predicted
# shape on day 7 and the red contour the predicted shape on day 150

contour <- vector('list')

for (i in 1:length(days)) {
  contour[[i]] <- exp_map(manifold, geodesics[[1]]$p, (days[i] -
      uncontam_mean_x) * geodesics[[1]]$V)
  contour[[i]] <- c(contour[[i]], contour[[i]][1])
}
plot(Re(contour[[length(days)]]), Im(contour[[length(days)]]), type = 'n',
    xaxt = 'n', yaxt = 'n', ann = FALSE, asp = 1)
 for (i in 1:length(days)) {
   lines(Re(contour[[i]]), Im(contour[[i]]), col = pal[i])
}
for (j in 2:4) {
  for (i in 1:length(days)) {
    contour[[i]] <- exp_map(manifold, geodesics[[j]]$p, (days[i] -
        contam_mean_x) * geodesics[[j]]$V)
    contour[[i]] <- c(contour[[i]], contour[[i]][1])
  }
  plot(Re(contour[[length(days)]]), Im(contour[[length(days)]]), type = 'n',
      xaxt = 'n', yaxt = 'n', ann = FALSE, asp = 1)
  for (i in 1:length(days)) {
    lines(Re(contour[[i]]), Im(contour[[i]]), col = pal[i])
  }
}
# even with a mere 4 corrupted landmarks out of a total of 8 * 168 = 1344, we
# can clearly see that contam_l2, the second image, looks slightly
# different from all the others, especially near the top of the image.

par(oldpar)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("exp_map")
### * exp_map

flush(stderr()); flush(stdout())

### Name: exp_map
### Title: Exponential map
### Aliases: exp_map

### ** Examples

exp_map('sphere', c(1, 0, 0, 0, 0), c(0, 0, pi / 4, 0, 0))




cleanEx()
nameEx("geo_dist")
### * geo_dist

flush(stderr()); flush(stdout())

### Name: geo_dist
### Title: Geodesic distance between two points on a manifold
### Aliases: geo_dist

### ** Examples

p1 <- matrix(rnorm(10), ncol = 2)
p1 <- p1[, 1] + (1i) * p1[, 2]
p1 <- (p1 - mean(p1)) / norm(p1 - mean(p1), type = '2')
p2 <- matrix(rnorm(10), ncol = 2)
p2 <- p2[, 1] + (1i) * p2[, 2]
p2 <- (p2 - mean(p2)) / norm(p2 - mean(p2), type = '2')
geo_dist('kendall', p1, p2)




cleanEx()
nameEx("geo_reg")
### * geo_reg

flush(stderr()); flush(stdout())

### Name: geo_reg
### Title: Gradient descent for (robust) geodesic regression
### Aliases: geo_reg

### ** Examples

# an example of multiple regression with two independent variables, with 64
# data points

x <- matrix(runif(2 * 64), ncol = 2)
x <- t(t(x) - colMeans(x))
y <- matrix(0L, nrow = 4, ncol = 64)
for (i in 1:64) {
  y[, i] <- exp_map('sphere', c(1, 0, 0, 0), c(0, runif(1), runif(1),
      runif(1)))
}
geo_reg('sphere', x, y, 'tukey', c = are_nr('tukey', 2, 6))




cleanEx()
nameEx("intrinsic_location")
### * intrinsic_location

flush(stderr()); flush(stdout())

### Name: intrinsic_location
### Title: Gradient descent for location based on M-type estimators
### Aliases: intrinsic_location

### ** Examples

y <- matrix(runif(100, 1000, 2000), nrow = 10)
intrinsic_location('euclidean', y, 'l2')




cleanEx()
nameEx("log_map")
### * log_map

flush(stderr()); flush(stdout())

### Name: log_map
### Title: Logarithm map
### Aliases: log_map

### ** Examples

log_map('sphere', c(0, 1, 0, 0), c(0, 0, 1, 0))




cleanEx()
nameEx("loss")
### * loss

flush(stderr()); flush(stdout())

### Name: loss
### Title: Loss
### Aliases: loss

### ** Examples

y <- matrix(0L, nrow = 3, ncol = 64)
for (i in 1:64) {
  y[, i] <- exp_map('sphere', c(1, 0, 0), c(0, runif(1), runif(1)))
}
intrinsic_mean <- intrinsic_location('sphere', y, 'l2')
loss('sphere', intrinsic_mean, numeric(3), numeric(64), y, 'l2')




cleanEx()
nameEx("onmanifold")
### * onmanifold

flush(stderr()); flush(stdout())

### Name: onmanifold
### Title: Manifold check and projection
### Aliases: onmanifold

### ** Examples

y1 <- matrix(rnorm(10), ncol = 2)
y1 <- y1[, 1] + (1i) * y1[, 2]
y2 <- matrix(rnorm(10), ncol = 2)
y2 <- y2[, 1] + (1i) * y2[, 2]
y3 <- matrix(rnorm(10), ncol = 2)
y3 <- y3[, 1] + (1i) * y3[, 2]
y3 <- (y3 - mean(y3)) / norm(y3 - mean(y3), type = '2') # project onto preshape space
y <- matrix(c(y1, y2, y3), ncol = 3)
onmanifold('kendall', y)




cleanEx()
nameEx("par_trans")
### * par_trans

flush(stderr()); flush(stdout())

### Name: par_trans
### Title: Parallel transport
### Aliases: par_trans

### ** Examples

p1 <- matrix(rnorm(10), ncol = 2)
p1 <- p1[, 1] + (1i) * p1[, 2]
p1 <- (p1 - mean(p1)) / norm(p1 - mean(p1), type = '2') # project onto pre-shape space
p2 <- matrix(rnorm(10), ncol = 2)
p2 <- p2[, 1] + (1i) * p2[, 2]
p2 <- (p2 - mean(p2)) / norm(p2 - mean(p2), type = '2') # project onto pre-shape space
p3 <- matrix(rnorm(10), ncol = 2)
p3 <- p3[, 1] + (1i) * p3[, 2]
p3 <- (p3 - mean(p3)) / norm(p3 - mean(p3), type = '2') # project onto pre-shape space
v <- log_map('kendall', p1, p3)
par_trans('kendall', p1, p2, v)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
