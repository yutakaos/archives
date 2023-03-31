#------------------------------------------------------------------------------#
# Script for main text
#------------------------------------------------------------------------------#
# Set working directory
setwd("2023_uic")
if (!dir.exists("figs")) dir.create("figs")

# Load library
library(rEDM)
library(rUIC)
library(RTransferEntropy)
library(pROC)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Read functions
source('functions/make_data.R')
source('functions/ccm_wrap.R')

#------------------------------------------------------------------------------#
# Figure S1
# Test function
test = function (block, lib="Y", tar="X", tp=1)
{
    out <- rUIC::uic.optimal(block, lib_var=lib, tar_var=tar, E=1:5, tp=-tp)
    E <- out$E[1]  # Optimal E for CCM and UIC
    data.frame(
        TE  = calc_te(block[,tar], block[,lib], lx=tp, ly=tp, shuffles=1),
        CCM = ccm_wrap(block, lib_var=lib, tar_var=tar, E=E, tp=-tp)$te[1],
        UIC = out$te[1]
    )
}

# AUC function for nonlinear dynamical system
get_auc_1 = function (tl = 200, zLib = 0.1, zTar = 0.1, Nrand = 1000)
{
    # data generations
    out_1 <- lapply(1:Nrand, function(k) {
        B <- make_data_2sp_logistic(tl=tl, Bxy=0.1, Zx=zTar, Zy=zLib)
        data.frame(causal=1, test(B))
    })
    out_0 <- lapply(1:Nrand, function(k) {
        B <- make_data_2sp_logistic(tl=tl, Bxy=0.0, Zx=zTar, Zy=zLib)
        data.frame(causal=0, test(B))
    })
    out <- do.call(rbind, c(out_1, out_0))
    out[out < 0] <- 0
    
    # compute AUC
    method <- colnames(out)[-1]
    roc <- lapply(method, function(x) {
        out <- pROC::roc(out[,1], out[,x], direction="<", quiet=TRUE)
        out <- with(out, data.frame(method=x, FPR=1-specificities, TPR=sensitivities))
        out <- out[nrow(out):1,]
        out
    })
    roc <- do.call(rbind, roc)
    roc$method <- factor(roc$method, levels=method)
    data.frame(
        tl=tl, zLib=zLib, zTar=zTar, roc,
        lib=sprintf("Effect's noise = %s", zLib),
        tar=sprintf("Cause's noise = %s" , zTar)
    )
}

# AUC function for linear dynamical system
get_auc_2 = function (tl = 200, zLib = 0.1, zTar = 0.1, Nrand = 1000)
{
    # data generations
    out_1 <- lapply(1:Nrand, function(k) {
        B <- make_data_2sp_var(tl=tl, Bxy=0.5, Ex=zTar, Ey=zLib)
        data.frame(causal=1, test(B))
    })
    out_0 <- lapply(1:Nrand, function(k) {
        B <- make_data_2sp_var(tl=tl, Bxy=0.0, Ex=zTar, Ey=zLib)
        data.frame(causal=0, test(B))
    })
    out <- do.call(rbind, c(out_1, out_0))
    out[out < 0] <- 0
    
    # compute AUC
    method <- colnames(out)[-1]
    roc <- lapply(method, function(x) {
        out <- pROC::roc(out[,1], out[,x], direction="<", quiet=TRUE)
        out <- with(out, data.frame(method=x, FPR=1-specificities, TPR=sensitivities))
        out <- out[nrow(out):1,]
        out
    })
    roc <- do.call(rbind, roc)
    roc$method <- factor(roc$method, levels=method)
    data.frame(
        tl=tl, zLib=zLib, zTar=zTar, roc,
        lib=sprintf("Effect's noise = %s", zLib),
        tar=sprintf("Cause's noise = %s" , zTar)
    )
}

# AUC function for nonlinear dynamical system with system noise
get_auc_3 = function (tl = 200, zLib = 0.1, zTar = 0.1, Nrand = 1000)
{
    # data generations
    out_1 <- lapply(1:Nrand, function(k) {
        B <- make_data_2sp_richer(tl=tl, Bxy=0.1, Ex=zTar, Ey=zLib)
        data.frame(causal=1, test(B))
    })
    out_0 <- lapply(1:Nrand, function(k) {
        B <- make_data_2sp_richer(tl=tl, Bxy=0.0, Ex=zTar, Ey=zLib)
        data.frame(causal=0, test(B))
    })
    out <- do.call(rbind, c(out_1, out_0))
    out[out < 0] <- 0
    
    # compute AUC
    method <- colnames(out)[-1]
    roc <- lapply(method, function(x) {
        out <- pROC::roc(out[,1], out[,x], direction="<", quiet=TRUE)
        out <- with(out, data.frame(method=x, FPR=1-specificities, TPR=sensitivities))
        out <- out[nrow(out):1,]
        out
    })
    roc <- do.call(rbind, roc)
    roc$method <- factor(roc$method, levels=method)
    data.frame(
        tl=tl, zLib=zLib, zTar=zTar, roc,
        lib=sprintf("Effect's noise = %s", zLib),
        tar=sprintf("Cause's noise = %s" , zTar)
    )
}

# Numerical simulations
# Warning messages are recieved from CCM, which has no effects on our results.
set.seed(870204)
df_1 <- rbind(
    get_auc_1(tl=100, zLib=0.01, zTar=0.01, Nrand=1000),
    get_auc_1(tl=100, zLib=0.01, zTar=0.5 , Nrand=1000),
    get_auc_1(tl=100, zLib=0.5 , zTar=0.01, Nrand=1000),
    get_auc_1(tl=100, zLib=0.5 , zTar=0.5 , Nrand=1000)
)
df_2 <- rbind(
    get_auc_2(tl=100, zLib=1, zTar=1, Nrand=1000),
    get_auc_2(tl=100, zLib=1, zTar=2, Nrand=1000),
    get_auc_2(tl=100, zLib=2, zTar=1, Nrand=1000),
    get_auc_2(tl=100, zLib=2, zTar=2, Nrand=1000)
)
df_3 <- rbind(
    get_auc_2(tl=100, zLib=0.1, zTar=0.1, Nrand=1000),
    get_auc_2(tl=100, zLib=0.1, zTar=1.0, Nrand=1000),
    get_auc_2(tl=100, zLib=1.0, zTar=0.1, Nrand=1000),
    get_auc_2(tl=100, zLib=1.0, zTar=1.0, Nrand=1000)
)
save(df_1, df_2, df_3, file="04_outs.RData")

# ggplot
gp_1 <- ggplot(df_1, aes(x=FPR, y=TPR, color=method)) +
    geom_line() + geom_abline(intercept=0, slope=1) +
    facet_grid(lib ~ tar) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    labs(x="FPR", y="TPR") +
    theme_bw(); print(gp_1)

gp_2 <- ggplot(df_2, aes(x=FPR, y=TPR, color=method)) +
    geom_line() + geom_abline(intercept=0, slope=1) +
    facet_grid(lib ~ tar) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    labs(x="FPR", y="TPR") +
    theme_bw(); print(gp_2)

gp_3 <- ggplot(df_3, aes(x=FPR, y=TPR, color=method)) +
    geom_line() + geom_abline(intercept=0, slope=1) +
    facet_grid(lib ~ tar) +
    coord_cartesian(xlim=c(0,1), ylim=c(0,1)) +
    labs(x="FPR", y="TPR") +
    theme_bw(); print(gp_3)

ggsave("figs/Fig_S1a.png", plot=gp_1, width=4.8, height=4.0)
ggsave("figs/Fig_S1b.png", plot=gp_2, width=4.8, height=4.0)
ggsave("figs/Fig_S1c.png", plot=gp_3, width=4.8, height=4.0)

# End