#------------------------------------------------------------------------------#
# 2-species nonlinear logistic map
#------------------------------------------------------------------------------#

make_data_2sp_logistic = function (tl = 400, Bxy = 0.1, Zx = 0, Zy = 0)
{
    x <- y <- rep(NA, tl)
    x[1] <- runif(1, 0.1, 0.5)
    y[1] <- runif(1, 0.1, 0.5)
    for (t in 1:(tl-1)) {
        x[t+1] <- x[t] * (3.8 - 3.8 * x[t] - 0.02 * y[t])
        y[t+1] <- y[t] * (3.5 - 3.5 * y[t] - Bxy  * x[t])
    }
    op <- data.frame(t=1:tl, x=x, y=y)
    op$X <- op$x + Zx * rnorm(tl, 0, sd(op$x, na.rm=TRUE))
    op$Y <- op$y + Zy * rnorm(tl, 0, sd(op$y, na.rm=TRUE))
    op
}

#------------------------------------------------------------------------------#
# 2-species linear VAR model
#------------------------------------------------------------------------------#

make_data_2sp_var = function (tl = 400, Bxy = 0.5, Ex = 2, Ey = 2)
{
    x <- y <- rep(0, tl + 1)
    for (t in seq(tl)) {
        x[t + 1] <- 0.1 * x[t]              + Ex * rnorm(1, 0, 1)
        y[t + 1] <- Bxy * x[t] + 0.1 * y[t] + Ey * rnorm(1, 0, 1)
    }
    data.frame(t=1:tl, X=x[-1], Y=y[-1])
}

#----------------------------------------------------------------------------------------------------#
# 2-species nonlinear Richer map with system noise
#----------------------------------------------------------------------------------------------------#

make_data_2sp_richer = function (tl = 400, Bxy = 0.1, Ex = 0, Ey = 0, Z = 0.0)
{
    x <- y <- rep(NA, tl)
    x[1] <- runif(1, 0.1, 0.5)
    y[1] <- runif(1, 0.1, 0.5)
    for (t in 1:(tl-1)) {
        x[t+1] <- x[t] * exp(3.8 - 3.8 * x[t] - 0.02 * y[t] + Ex * rnorm(1))
        y[t+1] <- y[t] * exp(3.5 - 3.5 * y[t] - Bxy  * x[t] + Ey * rnorm(1))
    }
    op <- data.frame(t = 1:tl, x = x, y = y)
    op$X <- op$x + Z * rnorm(tl, 0, sd(op$x, na.rm = TRUE))
    op$Y <- op$y + Z * rnorm(tl, 0, sd(op$y, na.rm = TRUE))
    op
}

# End