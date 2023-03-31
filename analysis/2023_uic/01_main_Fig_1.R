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
# Figure 1
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
    require(pROC)
    method <- colnames(out)[-1]
    roc <- lapply(method, function(x) roc(out[,1], out[,x], direction="<", quiet=TRUE))
    data.frame(
        tl=tl, zLib=zLib, zTar=zTar, method=factor(method, levels=method),
        lib=sprintf("Noise in X = %s", zLib), # Effect
        tar=sprintf("Noise in Y = %s", zTar), # Cause
        auc=sapply(roc, auc)
    )
}

# AUC function for linear dynamical system
get_auc_2 = function (tl = 200, zLib = 1, zTar = 1, Nrand = 1000)
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
    require(pROC)
    method <- colnames(out)[-1]
    roc <- lapply(method, function(x) roc(out[,1], out[,x], direction="<", quiet=TRUE))
    data.frame(
        tl=tl, zLib=zLib, zTar=zTar, method=factor(method, levels=method),
        lib=sprintf("Noise in X = %s", zLib), # Effect
        tar=sprintf("Noise in Y = %s", zTar), # Cause
        auc=sapply(roc, auc)
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
df_3 <- do.call(rbind, lapply(seq(10,300,10), get_auc_1, Nrand=1000))
df_4 <- do.call(rbind, lapply(seq(10,300,10), get_auc_2, Nrand=1000))
#save(df_1, df_2, df_3, df_4, file="01_outs.RData")

# ggplot
gp_1 <- ggplot(df_1, aes(x=method, y=auc, fill=method)) +
    geom_bar(stat="identity") +
    facet_grid(lib ~ tar) +
    coord_cartesian(ylim=c(0.5,1)) +
    labs(x=NULL, y="AUC") +
    theme_bw(); print(gp_1)

gp_2 <- ggplot(df_2, aes(x=method, y=auc, fill=method)) +
    geom_bar(stat="identity") +
    facet_grid(lib ~ tar) +
    coord_cartesian(ylim=c(0.5,1)) +
    labs(x=NULL, y="AUC") +
    theme_bw(); print(gp_2)

gp_3 <- ggplot(df_3, aes(x=tl, y=auc, color=method)) +
    geom_point() +
    stat_smooth(method="gam", formula=y ~ s(x, bs="cs"), se=FALSE) +
    coord_cartesian(ylim = c(0.5,1)) +
    labs(x="Number of Time Points", y="AUC") +
    theme_bw(); print(gp_3)

gp_4 <- ggplot(df_4, aes(x=tl, y=auc, color=method)) +
    geom_point() +
    stat_smooth(method="gam", formula=y ~ s(x, bs="cs"), se=FALSE) +
    coord_cartesian(ylim = c(0.5,1)) +
    labs(x="Number of Time Points", y="AUC") +
    theme_bw(); print(gp_4)

ggsave("figs/Fig_1a.png", plot=gp_1, width=4.8, height=4.0)
ggsave("figs/Fig_1b.png", plot=gp_2, width=4.8, height=4.0)
ggsave("figs/Fig_1c.png", plot=gp_3, width=4.8, height=2.4)
ggsave("figs/Fig_1d.png", plot=gp_4, width=4.8, height=2.4)

# End