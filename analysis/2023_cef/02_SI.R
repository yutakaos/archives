#------------------------------------------------------------------------------#
# Script for supplementary information materials
#------------------------------------------------------------------------------#
# Set working directory
setwd("2023_cef")
if (!dir.exists("figs")) dir.create("figs")

# Load library
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Read functions
source('functions/pred_EF.R')
source('functions/sims_EF.R')

#------------------------------------------------------------------------------#
# Figure S1
make_df = function (R = 10000)
{
    N = 10; S = 10*1:5; C = c(0.1,0.2,0.3)
    p = c(0.25,0.5,0.75)
    PX = c("Competitive", "Mixed", "Mutualistic")
    EF = c("Mean", "BC variability")
    
    compute_EF = function (x)
    {
        S <- length(x[[1]])
        ef_N <- sapply(x, function(x) sum(x[1:N]))
        ef_S <- sapply(x, function(x) sum(x))
        sd_x <- sapply(x, function(x) sd(x))
        fc   <- sapply(x, function(x) all(x>0))
        
        mEF1 <- mean(ef_N); vEF1 = sd(ef_N)/mEF1
        mEF2 <- mean(ef_S); vEF2 = sd(ef_S)/mEF2
        ef_o <- c(mEF1, vEF1, mEF2, vEF2)
        mEF1 <- mean(ef_N[fc]); vEF1 = sd(ef_N[fc])/mEF1
        mEF2 <- mean(ef_S[fc]); vEF2 = sd(ef_S[fc])/mEF2
        ef_x <- c(mEF1, vEF1, mEF2, vEF2)
        data.frame(ef1 = ef_o, ef2 = ef_x)
    }
    
    # Simulate CEF
    pars <- expand.grid(S=S, c=C, p=p)
    df_sims <- lapply(1:nrow(pars), function (i) {
        B <- with(pars[i,], lapply(1:R, sims_EF, S=S, c=c, p=p, apr=TRUE))
        x <- lapply(B, function(B) B$x)
        y <- lapply(B, function(B) B$y)
        out <- cbind(compute_EF(x), compute_EF(y))
        colnames(out) <- c("val","valF","valA","valAF")
        with(pars[i,], cbind(
            S=S, c=c, p=p, N=c("N","N","S","S"), out,
            ef=rep(EF,2), px=PX[4*p]))
    })
    df_sims <- do.call(rbind, df_sims)
    
    # Predict CEF
    pars <- expand.grid(S=min(S):max(S), c=C, p=p)
    df_pred <- lapply(1:nrow(pars), function(i) {
        S <- pars$S[i]; c = pars$c[i]; p = pars$p[i]
        val <- unlist(c(pred_EF(S,c,p,N), pred_EF(S,c,p,S)))
        names(val) <- NULL
        data.frame(
            S=S, c=c, p=p, N=c("N","N","S","S"), val=val,
            ef=rep(EF,2), px=PX[4*p])
    })
    df_pred <- do.call(rbind, df_pred)
    
    df_sims$ef <- factor(df_sims$ef, levels=EF)
    df_sims$px <- factor(df_sims$px, levels=PX)
    df_sims$c  <- factor(df_sims$c)
    df_pred$ef <- factor(df_pred$ef, levels=EF)
    df_pred$px <- factor(df_pred$px, levels=PX)
    df_pred$c  <- factor(df_pred$c)
    list(
        sims_N=df_sims[df_sims$N=="N",],
        pred_N=df_pred[df_pred$N=="N",],
        sims_S=df_sims[df_sims$N=="S",],
        pred_S=df_pred[df_pred$N=="S",]
    )
}
df = make_df(R = 1e6)
#save.image("out_SI.RData")

gp_1 <- ggplot(df$pred_S, aes(x = S, y = val, color = c)) +
    geom_line() +
    geom_line (data=df$sims_S, aes(y=valA), linetype="dashed") +
    geom_point(data=df$sims_S, aes(y=val)) +
    geom_point(data=df$sims_S, aes(y=valF), shape = 8) +
    facet_grid(ef ~ px, scales="free") +
    scale_color_discrete(name="Connectance", guide="none") +
    labs(x="Number of species", y="Ecosystem functioning") +
    theme_bw(); gp_1

gp_2 <- ggplot(df$pred_N, aes(x = S, y = val, color = c)) +
    geom_line() +
    geom_line (data=df$sims_N, aes(y=valA), linetype="dashed") +
    geom_point(data=df$sims_N, aes(y=val)) +
    geom_point(data=df$sims_N, aes(y=valF), shape = 8) +
    facet_grid(ef ~ px, scales="free") +
    scale_color_discrete(name="Connectance", guide="none") +
    labs(x="Number of species", y="Ecosystem functioning") +
    theme_bw(); gp_2
    
ggsave("figs/Fig_S1a.png", plot=gp_1, width=4.8, height=3.2)
ggsave("figs/Fig_S1b.png", plot=gp_2, width=4.8, height=3.2)

#------------------------------------------------------------------------------#
# Figure S2
make_df = function (S = 50, C = 0.3, N = 10)
{
    # Predict CEF
    p  = c(0.25, 0.5, 0.75)
    PX = c("Competitive", "Mixed", "Mutualistic")
    EF = c("Mean", "BC variability")
    type = c("SA","SB","PR","PC")
    zsym = seq(-1, 1, by=0.01)
    pars <- expand.grid(p=p, z=zsym, Z=type)
    df <- lapply(1:nrow(pars), function(i) {
        fn  <- with(pars, sprintf("pred_EF(S,C,%s,N,%s=%s)", p[i], Z[i], z[i]))
        val <- unlist(eval(parse(text=fn)))
        names(val) <- NULL
        with(pars, data.frame(
            S=S, c=C, p=p[i], N=N, Z=Z[i], z=z[i], val=val, ef=EF, px=PX[4*p][i],
            inZ =paste(ifelse(Z[i] %in% c("SA","PC"), "Positive", "Negative"), "incoming"),
            outZ=paste(ifelse(Z[i] %in% c("SA","PR"), "Positive", "Negative"), "outgoing")
        ))
    })
    df <- do.call(rbind, df)
    df$ef <- factor(df$ef, levels=EF)
    df$px <- factor(df$px, levels=PX)
    df$inZ  <- factor(df$inZ)
    df$outZ <- factor(df$outZ)
    return(df)
}
df <- make_df(S=50, C=0.3, N=10)

gp_3 <- ggplot(df[df$ef=="Mean",], aes(x=z, y=val, color=px)) +
    geom_line() +
    facet_grid(outZ ~ inZ, scales="free") +
    scale_color_discrete(guide="none") +
    labs(x="Symmetry", y="Ecosystem functioning") +
    theme_bw(); gp_3

gp_4 <- ggplot(df[df$ef!="Mean",], aes(x=z, y=val, color=px)) +
    geom_line() +
    facet_grid(outZ ~ inZ, scales="free") +
    scale_color_discrete(guide="none") +
    labs(x="Symmetry", y="Ecosystem functioning") +
    theme_bw(); gp_4

ggsave("figs/Fig_S2a.png", plot=gp_3, width=3.2, height=3.2)
ggsave("figs/Fig_S2b.png", plot=gp_4, width=3.2, height=3.2)

# End