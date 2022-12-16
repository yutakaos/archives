#------------------------------------------------------------------------------#
# Script for main text
#------------------------------------------------------------------------------#
# Set working directory
setwd("2022_cef")
if (!dir.exists("figs")) dir.create("figs")

# Load library
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Read functions
source('functions/pred_EF.R')

#------------------------------------------------------------------------------#
# Figure 1
make_df = function ()
{
    # Predict C-EF
    N = 10; S = 10:50; C = c(0.1,0.2,0.3)
    p = c(0.25,0.5,0.75)
    PX = c("Competitive", "Mixed", "Mutualistic")
    EF = c("Productivity", "BC variability")
    pars <- expand.grid(S=S, c=C, p=p)
    df_pred <- lapply(1:nrow(pars), function(i) {
        S <- pars$S[i]
        c <- pars$c[i]
        p <- pars$p[i]
        val <- unlist(c(pred_EF(S,c,p,N), pred_EF(S,c,p,S)))
        names(val) <- NULL
        data.frame(
            S=S, c=c, p=p, N=c("N","N","S","S"), val=val,
            ef=rep(EF,2), px=PX[4*p])
    })
    df_pred <- do.call(rbind, df_pred)
    df_pred$ef <- factor(df_pred$ef, levels=EF)
    df_pred$px <- factor(df_pred$px, levels=PX)
    df_pred$c  <- factor(df_pred$c)
    list(N=df_pred[df_pred$N=="N",], S=df_pred[df_pred$N=="S",])
}
df <- make_df()

gp_1 <- ggplot(df$S) +
    geom_line (aes(x=S, y=val, color=c)) +
    facet_grid(ef ~ px, scales="free") +
    scale_color_discrete(name="Connectance", guide="none") +
    labs(x="Number of species", y="Ecosystem functioning") +
    theme_bw(); gp_1

gp_2 <- ggplot(df$N) +
    geom_line (aes(x=S, y=val, color=c)) +
    facet_grid(ef ~ px, scales="free") +
    scale_color_discrete(name="Connectance", guide="none") +
    labs(x="Number of species", y="Ecosystem functioning") +
    theme_bw(); gp_2

ggsave("figs/Fig_1a.png", plot=gp_1, width=4.8, height=3.2)
ggsave("figs/Fig_1b.png", plot=gp_2, width=4.8, height=3.2)

#------------------------------------------------------------------------------#
# Figure 2
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

ggsave("figs/Fig_2a.png", plot=gp_3, width=3.2, height=3.2)
ggsave("figs/Fig_2b.png", plot=gp_4, width=3.2, height=3.2)

# End