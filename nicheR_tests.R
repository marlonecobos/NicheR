
# Tester File -------------------------------------------------------------

# Bulid elliposed tester

# Test 1: basic run
range <- data.frame(min = c(10, 20),
                    max = c(20, 35))

ell1 <- build_ellipsoid(range = range,
                       range_level = 0.95,
                       level = 0.99,
                       cor_tilt = NULL,
                       n_points = 200)

print(ell1)

# Test 2: with a tilt

cor_tilt <- matrix(data = c( 1.0, -0.7,
                             -0.7, 1.0),
                   nrow = 2)


ell2 <- build_ellipsoid(range = range,
                        cor_tilt = cor_tilt)


print(ell2)


# ell2 <- ell1

# Grab boundary points
p1 <- ell1$boundary_points
p2 <- ell2$boundary_points


# Plot window
rx <- range(c(p1[,1], p2[,1], ell1$mu[1], ell1$mu[1]))
ry <- range(c(p1[,2], p2[,2], ell1$mu[2], ell1$mu[2]))


plot(NA, xlim = rx, ylim = ry, asp = 1,
     xlab = "x1", ylab = "x2",
     main = "Ellipsoids: untilted vs tilted")


# Draw outlines (your 2D boundary points are in order, so lines works)
lines(p1[,1], p1[,2], lwd = 2) # untilted
lines(p2[,1], p2[,2], lwd = 2, lty = 2) # tilted


# Add centroids
points(ell1$mu[1], ell1$mu[2], pch = 16, cex = 1.1)
points(ell1$mu[1], ell1$mu[2], pch = 4, cex = 1.2, lwd = 2)


legend("topright",
       legend = c("untilted (R = I)",
                  "tilted (cor_tilt != I)"),
       lty = c(1, 2), lwd = 2, bty = "n")



