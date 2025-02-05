#------------------------------------------------------------------------------#
# Script for main figure 3 and SI figure 1
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
library(patchwork)

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
# 1. UIC
set.seed(1234)
out_1 <- NULL
for (i in 1:nrow(taxa)) {
    out <- uic_direct(
        df, group="plot", lib_var=taxa$Taxa_ID[i], tar_var=taxa$Taxa_ID,
        cond_var=NULL, E=1:10, tp=-(1:7), alpha=0.05)
    out_1 <- rbind(out_1,
        data.frame(cause=match(out$target, taxa$Taxa_ID), effect=i, out))
    print(sprintf("i = %s", i))
}; rm(i, out)
rownames(out_1) = NULL

# 2. S-map
out_2a <- fun_smap(out_1[out_1$pval<0.05 ,])
out_2b <- fun_smap(out_1[out_1$seq_test>0,])
out_2a$stats$sign <- ifelse(out_2a$stats$smap>0, "Positive", "Negative")
out_2b$stats$sign <- ifelse(out_2b$stats$smap>0, "Positive", "Negative")
save.image("03_main_outs.RData")

gp_1a <- with(taxa, gg_mynet(out_2a$stats, Taxa_ID, size=freq))
gp_2a <- with(taxa, gg_mynet(out_2b$stats, Taxa_ID, size=freq))
gp_1b <- ggplot(out_2a$stability, aes(x=date, y=val, color=plot)) +
    geom_point() + scale_y_continuous(limits=c(1,3.5)) +
    labs(x="Date (month)", y="Stability") +
    theme_classic() + theme(legend.position="none")
gp_2b <- ggplot(out_2b$stability, aes(x=date, y=val, color=plot)) +
    geom_point() + scale_y_continuous(limits=c(1,3.5)) +
    labs(x="Date (month)", y="Stability") +
    theme_classic() + theme(legend.position="none")
gp_1 <- gp_1a/gp_1b + plot_layout(height=c(2,1)); print(gp_1)
gp_2 <- gp_2a/gp_2b + plot_layout(height=c(2,1)); print(gp_2)
ggsave("figs/Fig_3a_.png", plot=gp_1, width=4, height=6)
ggsave("figs/Fig_3b_.png", plot=gp_2, width=4, height=6)

df <- data.frame(
    test = c(rep("Unconditional test",2),rep("Conditional test",2)),
    sign = c("Negative","Positive","Negative","Positive"),
    count = c(table(out_2a$stats$sign),table(out_2b$stats$sign))
)
df$test <- factor(df$test, levels=c("Unconditional test","Conditional test"))
df$sign <- factor(df$sign, levels=c("Negative","Positive"))
gp_3 <- ggplot(df, aes(x=test, y=count, fill=sign, color=sign)) +
    geom_bar(stat="identity", position="dodge", alpha=0.5) +
    scale_y_continuous(expand=c(0,0)) +
    labs(x="", y="Number of interactions") +
    theme_classic() + theme(legend.position="none"); print(gp_3)
ggsave("figs/Fig_3c_.png", plot=gp_3, width=3, height=4)

#------------------------------------------------------------------------------#
# Time-series for ggplot
df <- data.frame(
    date = rep(smpl$date, nrow(taxa)),
    plot = rep(smpl$plot, nrow(taxa)),
    taxa = c(sapply(taxa$Taxa_ID, rep, length=nrow(smpl))),
    asv  = c(asv)
)
rownames(df) <- NULL

gp_3 <- ggplot(df, aes(x=date, y=asv, color=plot)) +
    geom_line() + facet_wrap(~ taxa) +
    labs(x="Date (month)", y="Standardized DNA copy") + 
    theme_bw() + theme(legend.position="none")
ggsave("figs/Fig_S1.png", plot=gp_3, width=7.2, height=7.2)

# End