#------------------------------------------------------------------------------#
# Script for main figure 4 and SI figure 1
#------------------------------------------------------------------------------#
# Set working directory
setwd("2026_uic")
if (!dir.exists("figs")) dir.create("figs")

# Load packages
library(rEDM)
library(rUIC)
library(RTransferEntropy)
library(foreach)
library(doParallel)
library(ggplot2)
library(patchwork)

# Read functions
source("functions/make_data_multisp.R")
source("functions/test_multisp.R")
source("functions/uic_direct.R")
rm(fun_smap, gg_mynet)

#------------------------------------------------------------------------------#
# Test function
test = function (
    ncl = 25, nsp = 40, model = c("Ricker","VAR"), N = 100)
{
    fun_test <<- function (data) {
        df <- data$df
        sp_name <- colnames(df[,-1])
        
        # TE
        out_te <- lapply(sp_name, function(x) te_test(df, lib_var=x, tar_var=sp_name))
        out_te <- do.call(rbind, out_te)
        M_te <- make_M_from_test(out_te, nsp=length(sp_name))
        
        # CCM
        out_ccm <- lapply(sp_name, function(x) ccm_test(df, lib_var=x, tar_var=sp_name, maxE=10))
        out_ccm <- do.call(rbind, out_ccm)
        M_ccm <- make_M_from_test(out_ccm, nsp=length(sp_name))
        
        # UIC
        out_uic <- lapply(sp_name, function(x) 
            uic_direct(df, lib_var=x, tar_var=sp_name, E=1:10, tp=-1))
        out_uic <- do.call(rbind, out_uic)
        out_uic$pval2 <- out_uic$pval
        out_uic$pval2[out_uic$seq_test<0] <- 1
        M_uic1 <- make_M_from_test(out_uic, nsp=length(sp_name))
        M_uic2 <- make_M_from_test(out_uic, nsp=length(sp_name), pval="pval2")
        
        # output
        out <- rbind(
            cbind(method="TE"  , calc_metrics(data$M, M_te)),
            cbind(method="CCM" , calc_metrics(data$M, M_ccm)),
            cbind(method="UIC" , calc_metrics(data$M, M_uic1)),
            cbind(method="cUIC", calc_metrics(data$M, M_uic2))
        )
        out$method <- factor(out$method, levels=c("TE","CCM","UIC","cUIC"))
        out
    }
    on.exit(rm(fun_test, envir=.GlobalEnv))
    
    # Set parallel computing
    model <- match.arg(model)
    cl <- makeCluster(ncl); registerDoParallel(cl)
    packages <- c("rEDM","rUIC","RTransferEntropy")
    export <- c(
        "make_data_multisp_ricker","make_data_multisp_var","fun_test",
        "ccm_test","te_test","uic_direct","calc_metrics","make_M_from_test")
    clusterExport(cl, export)
    on.exit(stopCluster(cl))
    
    # Test
    df <- data.frame()
    cat(sprintf("n = 0 : %s\n", Sys.time()))
    while(nrow(df)/4 < N) {
        kmax <- min(ncl, N-nrow(df)/4)
        df0 <- foreach(k=1:kmax+nrow(df)/4, .combine="rbind", .packages=packages) %dopar% {
            if (model=="Ricker") make_data = make_data_multisp_ricker
            if (model=="VAR")    make_data = make_data_multisp_var   
            data <- make_data(tl=100, nsp=nsp, connectance=0.2)
            cbind(no=k, model=model, nsp=nsp, fun_test(data))
        }
        df <- rbind(df, df0); rm(df0)
        cat(sprintf("n = %s : %s\n", nrow(df)/4, Sys.time()))
    }
    return(df)
}

#------------------------------------------------------------------------------#
# Test
set.seed(1234)
df_r20 <- test(ncl=25, nsp=20, model="Ricker", N=100)
df_v20 <- test(ncl=25, nsp=20, model="VAR"   , N=100)
df_r40 <- test(ncl=25, nsp=40, model="Ricker", N=100)
df_v40 <- test(ncl=25, nsp=40, model="VAR"   , N=100)
save.image("03_main_outs.RData")

# Plot 1
gg_fun = function(df, metric=c("shd","acc"))
{
    if (metric == "shd")
        gp <- ggplot(df) + 
            geom_boxplot(aes(x=method, y=shd, fill=method)) +
            labs(x=NULL, y="Structural Hamming distance")
    if (metric == "acc")
        gp <- ggplot(df) + 
            geom_boxplot(aes(x=method, y=acc, fill=method)) +
            labs(x=NULL, y="Accuracy")
    gp <- gp + theme_classic() + theme(legend.position="none")
    return(gp)
}
print(gp_1a <- gg_fun(df_r20, metric="acc") + ggtitle("(a) 20-species Richker"))
print(gp_1b <- gg_fun(df_r40, metric="acc") + ggtitle("(b) 40-species Richker"))
print(gp_1c <- gg_fun(df_v20, metric="acc") + ggtitle("(c) 20-species VAR"))
print(gp_1d <- gg_fun(df_v40, metric="acc") + ggtitle("(d) 40-species VAR"))
gp_1 <- gp_1a+gp_1b+gp_1c+gp_1d+plot_layout(nrow=2, ncol=2)
ggsave("figs/Fig_4.png", plot=gp_1, width=6, height=8)

print(gp_2a <- gg_fun(df_r20, metric="shd") + ggtitle("(a) 20-species Richker"))
print(gp_2b <- gg_fun(df_r40, metric="shd") + ggtitle("(b) 40-species Richker"))
print(gp_2c <- gg_fun(df_v20, metric="shd") + ggtitle("(c) 20-species VAR"))
print(gp_2d <- gg_fun(df_v40, metric="shd") + ggtitle("(d) 40-species VAR"))
gp_2 <- gp_2a+gp_2b+gp_2c+gp_2d+plot_layout(nrow=2, ncol=2)
ggsave("figs/Fig_S1.png", plot=gp_2, width=6, height=8)

# Plot 2
count_best_method = function (df)
{
    acc <- apply(matrix(df[,"acc"], nrow=4), 2, which.max)
    shd <- apply(matrix(df[,"shd"], nrow=4), 2, which.min)
    out <- data.frame(
        method = c("TE","CCM","UIC","cUIC"),
        acc = sapply(1:4, function(k) sum(acc==k)),
        shd = sapply(1:4, function(k) sum(shd==k))
    )
    out$method <- factor(out$method, levels=out$method)
    out
}
df <- rbind(
    cbind(data="20sp Ricker", count_best_method(df_r20)),
    cbind(data="40sp Ricker", count_best_method(df_r40)),
    cbind(data="20sp VAR"   , count_best_method(df_v20)),
    cbind(data="40sp VAR"   , count_best_method(df_v40))
)
df$data <- factor(df$data, levels=c("20sp Ricker","40sp Ricker","20sp VAR","40sp VAR"))
levels(df$data) <- c("20sp\nRicker","40sp\nRicker","20sp\nVAR","40sp\nVAR")
gp_1e <- ggplot(df, aes(x=data, y=acc, fill=method)) +
    geom_bar(stat="identity", position="fill") +
    scale_y_continuous(labels = scales::percent) +
    labs(x=NULL, y="Accucary", fill=NULL) +
    theme_bw() + ggtitle("(e) Proportion of best method"); gp_1e
ggsave("figs/Fig_4e.png", plot=gp_1e, width=3.5, height=4)

gp_2e <- ggplot(df, aes(x=data, y=shd, fill=method)) +
    geom_bar(stat="identity", position="fill") +
    scale_y_continuous(labels = scales::percent) +
    labs(x=NULL, y="Accucary", fill=NULL) +
    theme_bw() + ggtitle("(e) Proportion of best method"); gp_2e
ggsave("figs/Fig_S1e.png", plot=gp_2e, width=3.5, height=4)

# End