#------------------------------------------------------------------------------#
# Script for main text
#------------------------------------------------------------------------------#
# Set working directory
setwd("2023_uic")
if (!dir.exists("figs")) dir.create("figs")

# Load library
library(rUIC)
library(rEDM)
library(ggplot2)
library(igraph)
library(ggraph)
library(cowplot)
theme_set(theme_cowplot())

# Read functions
source('functions/uic_direct.R')

# Read data
asv  <- read.csv("ushio2022/data_asvtable.csv")
smpl <- read.csv("ushio2022/data_sample.csv")
taxa <- read.csv("ushio2022/data_taxa.csv")
smpl$date <- as.Date(smpl$date)
smpl$plot <- as.factor(smpl$plot)

# Extract only Proteobacteria with high occurence during July
taxa <- with(taxa, taxa[phylum=="Proteobacteria",])
asv  <- asv[,colnames(asv) %in% taxa$Taxa_ID]
freq <- apply(asv["2017-07-01"<=smpl$date&smpl$date<="2017-07-31",]>0, 2, mean)
taxa$freq <- freq; rm(freq)
taxa <- with(taxa, taxa[freq>0.5,])
asv  <- asv[,colnames(asv) %in% taxa$Taxa_ID]
asv  <- apply(asv, 2, function(x) x/max(x))
df   <- data.frame(date=smpl$date, plot=smpl$plot, asv)

#------------------------------------------------------------------------------#
# 1. CCM
out_1 <- NULL
for (i in 1:nrow(taxa)) {
    out <- uic_direct(
        df, group="plot", lib_var=taxa$Taxa_ID[i], tar_var=taxa$Taxa_ID,
        cond_var=NULL, E=1:10, tp=-(1:7))
    out_1 <- rbind(out_1,
        data.frame(cause=match(out$target, taxa$Taxa_ID), effect=i, out))
    print(sprintf("i = %s", i))
}; rm(i, out)
rownames(out_1) = NULL

# 2. S-map
fun_smap = function (out)
{
    # Make delay time-series
    make_dc = function (X, lib = c(1,NROW(X)), E = 0, group=NULL)
    {
        lib <- rbind(lib)
        if (is.null(group)) {
            X <- as.matrix(X)
            NAsF <- matrix(NA, nrow=pmax(0, E), ncol=ncol(X))
            NAsB <- matrix(NA, nrow=pmax(0,-E), ncol=ncol(X))
            nNA <- nrow(NAsB)
            out <- lapply(1:nrow(lib), function(i) {
                L1 <- lib[i,1]
                L2 <- lib[i,2]
                dc <- rbind(NAsF, X[L1:L2,,drop=FALSE], NAsB)
                dc [nNA+1:(L2-L1+1),,drop=FALSE]
            })
            out <- do.call(rbind, out)
            return(out)
        }
        lib <- lapply(1:nrow(lib), function(i) {
            L <- lapply(unique(group), function(g) range((lib[i,1]:lib[i,2])[group==g]))
            L <- do.call(rbind, L)
        })
        lib <- do.call(rbind, lib)
        make_dc(X, lib, E)
    }
    
    smap <- lapply(1:nrow(taxa), function(i) {
        # theta to be used for checking nonlinearty
        theta <- c(
            0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5,
            0.75, 1, 1.5, 2, 3, 4, 6, 8)
        
        # make block
        outi <- out[out$effect==i,]
        lib_var <- taxa$Taxa_ID[i]
        if (nrow(outi) != 0) {
            B <- lapply(1:nrow(outi), function(k) {
                x <- taxa$Taxa_ID[outi$cause[k]]
                make_dc(df[,x], E=-outi$tp[k]-1, group=df$plot)
            })
            B <- data.frame(do.call(cbind, B))
            cond_var <- with(outi, paste0(taxa$Taxa_ID[cause],".",-tp-1))
            colnames(B) <- cond_var
        }
        else {
            B <- matrix(NA, nrow=nrow(df), ncol=0)
            cond_var <- NULL
        }
        
        # simplex
        simp <- rUIC::simplex(
            data.frame(plot=df$plot, effect=df[,lib_var], B), group="plot",
            lib_var="effect", cond_var=cond_var, E=1:10, tp=1, alpha=0.05)
        optE <- with(simp, max(0, E[pval < 0.05]))
        
        # smap
        lib <- lapply(unique(df$plot), function(x) range(which(df$plot==x)))
        lib <- do.call(rbind, lib)
        L <- lapply(0:optE-1, function(k) make_dc(df[,lib_var], E=k, group=df$plot))
        L <- data.frame(do.call(cbind, L))
        n <- ncol(L)+ncol(B)
        smap <- rEDM::block_lnlp(
            cbind(L,B), lib=lib, method="s-map", tp=1,
            columns= 2:n, target_column=1, theta=theta,
            silent=TRUE, save_smap_coefficients=TRUE)
        coef <- smap$smap_coefficients[which.min(smap$rmse)][[1]]
        coef <- cbind(coef[,n,drop=FALSE], coef[,-n,drop=FALSE])
        med  <- apply(coef[,-1:-ncol(L),drop=FALSE], 2, median, na.rm=TRUE)
        smry <- data.frame(
            cause  = c(rep(lib_var, optE+1), taxa$Taxa_ID[outi$cause]),
            effect = lib_var, E = c(0:optE, -outi$tp))
        
        # output
        list(stats = cbind(outi, smap=med), summary = smry, coef = coef)
    })
    
    # dynamical stability
    smry <- do.call(rbind, lapply(smap, function(x) x$summary))
    coef <- do.call(cbind, lapply(smap, function(x) x$coef))
    coef <- coef[,smry$E!=0]
    smry <- smry[smry$E!=0,]
    smry$i <- match(smry$cause , taxa$Taxa_ID)
    smry$j <- match(smry$effect, taxa$Taxa_ID)
    time <- which(complete.cases(coef))
    stab <- data.frame(df[,1:2], val=NA)
    stab$val[time] <- sapply(time, function(t) {
        nsp  <- nrow(taxa)
        adj  <- data.frame(smry, val=unlist(coef[t,]))
        Jt <- matrix(0, nrow=nsp, ncol=max(adj$E)*nsp)
        for (k in 1:nrow(adj)) {
            Jt[adj$i[k], (adj$E[k]-1)*nsp+adj$j[k]] <- adj$val[k]
        }; rm(k)
        Jt <- rbind(Jt, diag(nrow=ncol(Jt)-nsp, ncol=ncol(Jt)))
        abs(eigen(Jt)$values[1])
    })
    
    # output
    list(
        stats = do.call(rbind, lapply(smap, function(x) x$stats)),
        stability = stab
    )
}
out_2a <- fun_smap(out_1[out_1$pval<0.05  ,])
out_2b <- fun_smap(out_1[out_1$direct=="+",])
save.image("03_outs.RData")

# ggplot
gg_mynet = function (out, name, size=1)
{
    if (length(size)==1) size <- rep(size, length(name))
    mygraph <- igraph::graph_from_data_frame(
        data.frame(from="origin", to=name),
        vertices=data.frame(name=c("origin", name), value=c(NA, size)))
    from <- out$cause  + 1  # +1 (origin)
    to   <- out$effect + 1  # +1 (origin)
    ggraph(mygraph, layout='dendrogram', circular=TRUE) + 
        geom_node_point(aes(filter=leaf, size=value), color="coral") +
        geom_conn_bundle(
            data=get_con(from=from, to=to, col=factor(out$smap>0)),
            aes(x=x*0.95, y=y*0.95, color=col), alpha=0.5, tension=0.7,
            arrow=grid::arrow(angle=20, length=unit(0.1,"inches"), type="closed")) +
        theme_void() + theme(legend.position="none")
}

gp_1 <- with(taxa, gg_mynet(out_2a$stats, Taxa_ID, size=freq))
gp_2 <- with(taxa, gg_mynet(out_2b$stats, Taxa_ID, size=freq))
gp_3 <- ggplot(out_2a$stability, aes(x=date, y=val, color=plot)) +
    geom_point() + scale_y_continuous(limits=c(1,4.4)) +
    labs(x="Date (Month)", y="Stability") +
    theme_bw() + theme(legend.position="none")
gp_4 <- ggplot(out_2b$stability, aes(x=date, y=val, color=plot)) +
    geom_point() + scale_y_continuous(limits=c(1,4.4)) +
    labs(x="Date (Month)", y="Stability") +
    theme_bw() + theme(legend.position="none")

ggsave("figs/Fig_3a.png", plot=gp_1, width=3.6, height=3.6)
ggsave("figs/Fig_3b.png", plot=gp_2, width=3.6, height=3.6)
ggsave("figs/Fig_3c.png", plot=gp_3, width=3.6, height=2.4)
ggsave("figs/Fig_3d.png", plot=gp_4, width=3.6, height=2.4)

#------------------------------------------------------------------------------#
# Time-series for ggplot
df <- data.frame(
    date = rep(smpl$date, nrow(taxa)),
    plot = rep(smpl$plot, nrow(taxa)),
    taxa = c(sapply(taxa$Taxa_ID, rep, length=nrow(smpl))),
    asv  = c(asv)
)
rownames(df) <- NULL

gp_5 <- ggplot(df, aes(x=date, y=asv, color=plot)) +
    geom_line() + facet_wrap(~ taxa) +
    labs(x="Date (Month)", y="Standardized DNA copy") + 
    theme_bw() + theme(legend.position="none")

ggsave("figs/Fig_S2.png", plot=gp_5, width=7.2, height=7.2)

# End