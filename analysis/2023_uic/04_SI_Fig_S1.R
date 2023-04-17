#------------------------------------------------------------------------------#
# Script for main text
#------------------------------------------------------------------------------#
# Set working directory
setwd("2023_uic")
if (!dir.exists("figs")) dir.create("figs")

# Load library
library(rUIC)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# 4-species logistic model (Ye et al. 2015)
make_data = function (tl = 200, Z = 0, init_y = rep(0.4, 4))
{
    y <- matrix(NA, tl, 4)
    y[1,] <- init_y[1:4]
    for (t in 1:(tl-1)) {
        y[t+1,1] <- y[t,1] * (3.9 - 3.9 * y[t,1])
        y[t+1,2] <- y[t,2] * (3.6 - 3.6 * y[t,2] - 0.4  * y[t,1])
        y[t+1,3] <- y[t,3] * (3.6 - 3.6 * y[t,3] - 0.4  * y[t,2])
        y[t+1,4] <- y[t,4] * (3.8 - 3.8 * y[t,4] - 0.35 * y[t,3])
    }
    y[,1] <- y[,1] + Z * rnorm(tl, 0, sd(y[,1], na.rm = TRUE))
    y[,2] <- y[,2] + Z * rnorm(tl, 0, sd(y[,2], na.rm = TRUE))
    y[,3] <- y[,3] + Z * rnorm(tl, 0, sd(y[,3], na.rm = TRUE))
    y[,4] <- y[,4] + Z * rnorm(tl, 0, sd(y[,4], na.rm = TRUE))
    colnames(y) <- sprintf("y%s", 1:4)
    data.frame(t=1:tl, y)
}

#------------------------------------------------------------------------------#
#
# Simulations
set.seed(870204)
block <- make_data(tl=200, init_y=runif(4, 0.1, 0.5))
Y <- colnames(block)[-1]

gp <- rep(list(NULL), 9)
for (i in 1:3) for (j in 1:3) {
    if (i==1 && j==2) {
        block$cond <- c(NA,block[-200,"y2"])
        out <- rUIC::uic.optimal(
            block, lib_var="y3", tar_var="y1", cond_var="cond",
            alpha=0.01, E=1:10, tp=-(0:8))
    }
    else if (i==1 && j==3) {
        block$cond <- c(NA,block[-200,"y3"])
        out <- rUIC::uic.optimal(
            block, lib_var="y4", tar_var="y1", cond_var="cond",
            E=1:10, tp=-(0:8))
    }
    else if (i==2 && j==3) {
        block$cond <- c(NA,block[-200,"y3"])
        out <- rUIC::uic.optimal(
            block, lib_var="y4", tar_var="y2", cond_var="cond",
            E=1:10, tp=-(0:8))
    }
    else {
        out <- rUIC::uic.optimal(
            block, lib_var=Y[i+1], tar_var=Y[j],
            E=1:10, tp=-(0:8))
    }
    gp[[3*(j-1) + i]] <- ggplot(out, aes(x=-tp, y=te)) +
        geom_line() + geom_hline(yintercept=0) +
        geom_point(size=2, show.legend=FALSE) +
        ylim(c(-0.1, 1.7)) + labs(x=NULL, y=NULL) +
        theme_bw()
}; rm(i, j, out)

gp <- eval(parse(text=paste0(
    "plot_grid(", paste0(sprintf("gp[[%s]],", 1:9), collapse=""),
    "nrow=3, ncol=3)"))); print(gp)

ggsave("figs/Fig_S1.png", plot=gp, width=8, height=8)

# End