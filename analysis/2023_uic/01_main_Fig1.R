#------------------------------------------------------------------------------#
# Script for main figure 1
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
library(foreach)
library(doParallel)
library(patchwork)

# Read functions
source('functions/make_data.R')

# Set number of simulations
Nsim <- 1000

#------------------------------------------------------------------------------#
# Test function with 3 models
# 1. logistic model with observation noise
# 2. linear VAR model
# 3. Richer model with system noise
test = function (tl = 100, eLib = 0.1, eTar = eLib, N = 1000, model="")
{
    if (!model%in%c("Logistic","VAR","Richer")) stop("Invalid model.")
    if (model=="Logistic") { make_data <<- make_data_2sp_logistic; Bxy <- 0.1 }
    if (model=="VAR")      { make_data <<- make_data_2sp_var     ; Bxy <- 0.5 }
    if (model=="Richer")   { make_data <<- make_data_2sp_richer  ; Bxy <- 0.1 }
    calc_ce <<- function (block, lib="Y", tar="X", tp=1) {
        E <- EmbedDimension(
            dataFrame=block, lib=c(1,nrow(block)), pred=c(1,nrow(block)),
            maxE=5, Tp=tp, columns=lib, target=lib, showPlot=FALSE)
        E <- E$E[which.max(E$rho)]
        out_ccm <- CCM(
            dataFrame=block, E=E, Tp=-tp, columns=lib, target=tar,
            libSizes=c(E+5, nrow(block)-E+1), sample=100, includeData=TRUE)
            # E+5 is used as minimum library size to suppress warnings.
        out_ccm <- aggregate(RMSE~LibSize, data=out_ccm$CCM1_PredictStat, mean)
        out_uic <- uic.optimal(block, lib_var=lib, tar_var=tar, E=1:5, tp=-tp)
        data.frame(
            TE  = calc_te(block[,tar], block[,lib], lx=tp, ly=tp, shuffles=1),
            CCM = with(out_ccm, log(RMSE[1])-log(RMSE[2])),
            UIC = out_uic$te[1]
        )
    }
    on.exit(rm(make_data, calc_ce, envir=.GlobalEnv))
    packages <- c("rEDM","rUIC","RTransferEntropy")
    clusterExport(cl, c("make_data","calc_ce"))
    out_1 <- foreach(k=1:N, .packages=packages) %dopar% {
        B <- make_data(tl=tl, Bxy=Bxy, eTar, eLib)
        data.frame(causal=1, calc_ce(B))
    }
    out_0 <- foreach(k=1:N, .packages=packages) %dopar% {
        B <- make_data(tl=tl, Bxy=0.0, eTar, eLib)
        data.frame(causal=0, calc_ce(B))
    }
    out <- do.call(rbind, c(out_1,out_0))
    out[out<0] <- 0
    method <- colnames(out)[-1]
    roc <- lapply(method, function(x) roc(out[,1], out[,x], direction="<", quiet=TRUE))
    data.frame(
        tl=tl, eLib=eLib, eTar=eTar, method=factor(method, levels=method),
        lib=sprintf("Noise in X = %s", eLib), # Effect
        tar=sprintf("Noise in Y = %s", eTar), # Cause
        auc=sapply(roc, auc)
    )
}

# Numerical experiments
exp1 <- function(e0, e1, model)
{
    rbind(
        test(eLib=e0, eTar=e0, N=Nsim, model=model),
        test(eLib=e0, eTar=e1, N=Nsim, model=model),
        test(eLib=e1, eTar=e0, N=Nsim, model=model),
        test(eLib=e1, eTar=e1, N=Nsim, model=model)
    )
}
exp2 <- function (e0, model)
{
    out <- lapply(seq(20,300,10), test, eLib=e0, N=Nsim, model=model)
    do.call(rbind, out)
}

#------------------------------------------------------------------------------#
# Numerical simulations
set.seed(1234)
cl <- makeCluster(60); registerDoParallel(cl)
df_1a <- exp1(e0=0.05, e1=0.5, model="Logistic")
df_1b <- exp1(e0=1.0 , e1=2.0, model="VAR")
df_1c <- exp1(e0=0.05, e1=0.5, model="Richer")
df_2a <- exp2(e0=0.05, model="Logistic")
df_2b <- exp2(e0=1.0 , model="VAR")
df_2c <- exp2(e0=0.05, model="Richer")
stopCluster(cl); rm(cl)
save.image("01_main_outs.RData")

gg_fun <- function(df1, df2) {
    gp_1 <- ggplot(df1, aes(x=method, y=auc, fill=method)) +
        geom_bar(stat="identity") +
        facet_grid(lib ~ tar) +
        coord_cartesian(ylim=c(0.5,1)) +
        labs(x=NULL, y="AUC") +
        theme_bw() + theme(legend.position="none")
    gp_2 <- ggplot(df2, aes(x=tl, y=auc, color=method)) +
        geom_point() +
        stat_smooth(method="gam", formula=y ~ s(x, bs="cs", k=15), se=FALSE) +
        coord_cartesian(ylim = c(0.5,1)) +
        labs(x="Number of time points", y="AUC") +
        theme_classic() + theme(legend.position="none")
    out <- gp_1/gp_2 + plot_layout(height=c(2,1)); print(out)
    out
}
gp_1a <- gg_fun(df_1a, df_2a)
gp_1b <- gg_fun(df_1b, df_2b)
gp_1c <- gg_fun(df_1c, df_2c)
ggsave("figs/Fig_1a.png", plot=gp_1a, width=4, height=6.4)
ggsave("figs/Fig_1b.png", plot=gp_1b, width=4, height=6.4)
ggsave("figs/Fig_1c.png", plot=gp_1c, width=4, height=6.4)

# End