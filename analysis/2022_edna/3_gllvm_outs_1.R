#------------------------------------------------------------------------------#
# GLLVM analysis 1
#------------------------------------------------------------------------------#
# Set directory
setwd("D:/Desktop/2022_edna")
if (!dir.exists("fig")) dir.create("fig")

# Load library
library(pROC)
library(corrplot)
library(ggplot2)
library(cowplot)
library(maps)
theme_set(theme_cowplot())

# Read data
name <- read.csv("data/20211110_species.csv", header=TRUE)
site <- read.csv("data/20211110_site.csv"   , header=TRUE)
cm   <- read.csv("data/20211110_community_matrix.csv", header=TRUE)

# Model-input data
idxS <- name$input == 1
idxL <- site$input == 1
name <- name[idxS,]
site <- site[idxL,]
cm <- cm[idxS, idxL]
rm(idxS, idxL)

# Map Polygon
Jpn <- map_data("world", region="Japan")

# Load RData
load("out/op_lvm_f3_g2_ests.RData")
ne = ncol(ests$mean$beta) - 1
nL = nrow(ests$mean$scoreR)
nS = nrow(ests$mean$psiR)

#------------------------------------------------------------------------------#
# In-sample AUC comparison
compute_AUC = function (nf = 0, ng = 0)
{
    # Load RData
    load(sprintf("out/op_lvm_f%s_g%s_ests.RData", nf, ng))
    expl  <- apply(site[,c("temperature","salinity")], 2, scale)
    coef  <- with(ests$mean, cbind(beta, psiR, psiS, 1))
    score <- with(ests$mean, cbind(1, expl, scoreR, scoreS, alpha))
    # Compute AUC
    mu_all <- coef %*% t(as.matrix(score))
    cm_all <- as.matrix(1 * (cm > 0))
    dimnames(mu_all) <- dimnames(cm_all) <- NULL
    sapply(1:nrow(mu_all), function(i)
        auc(roc(cm_all[i,] ~ pnorm(mu_all[i,]), quiet=TRUE))
    )
}
df <- data.frame(
    best = compute_AUC(nf=3, ng=2),
    null = compute_AUC(nf=0, ng=0),
    n = rowSums(cm > 0)
)

gp_1 <- ggplot(df, aes(x=null, y=best, color=n)) +
    geom_point() +
    scale_color_gradient(low="deepskyblue", high="deeppink") +
    geom_abline(intercept=0, slope=1) +
    xlim(c(0.5, 1.0)) + ylim(c(0.5, 1.0)) +
    labs(color = "Occurrence") + theme_bw() +
    xlab("In-sample AUC of null model") + ylab("In-sample AUC of best model")

ggsave("fig/Fig_S3_AUC.png", gp_1, width=6, height=6)
rm(compute_AUC, df)

#------------------------------------------------------------------------------#
# Map latent variables
lv <- cbind(ests$mean$scoreR, ests$mean$scoreS)
df <- lapply(1:(ng+nf), function(k) {
    clv = ifelse(k <= ng, "Local LV", "Regional LV")
    klv = ifelse(k <= ng, k, k - ng)
    data.frame(site[,c("lat","lon")], val=lv[,k], clv=clv, klv=klv)
})
df <- do.call(rbind, df)

# Explained variance
coef <- with(ests$mean, cbind(beta[,-1], psiR, psiS))
var <- apply(coef^2, 2, sum) / (nS + 1 - c(1,1,1:ng,1:nf))
var <- round(100 * var / sum(var), 1)

# Range parameters of Matern distance
rho <- round(sqrt(8) / ests$mean$kappa, 1)

# Annotations
annots <- data.frame(
    clv = c(rep("Local LV", ng), rep("Regional LV", nf)),
    klv = c(1:ng, 1:nf),
    stats = c(
        paste0("Explained = ", var[ne+1:ng], "%"),
        paste0("Explained = ", var[ne+ng+1:nf], "%\nRange = ", rho, "km"))
)

gp_2 <- ggplot(df) +
    coord_fixed(xlim=range(Jpn$long), ylim=range(Jpn$lat)) +
    geom_polygon(aes(x=long, y=lat, group=group, fill=region), data=Jpn, fill="grey70") +
    geom_point(aes(x=lon, y=lat, color=val), size=0.9) +
    scale_color_gradient2(low="deeppink", mid="beige", high="deepskyblue", midpoint=0) +
    scale_y_continuous(position="right") +
    facet_grid(clv ~ klv, switch="y") +
    geom_text(data=annots, aes(x=124, y=45, label=stats), size=3, hjust="left", vjust="top") +
    labs(color="Score") +
    theme_bw() + xlab("Longitude") + ylab("Latitude")

ggsave("fig/Fig_2_LV_Map.png", gp_2, width=9, height=6)
rm(lv, df, coef, var, rho, annots)

#------------------------------------------------------------------------------#
# Correlation plot
calc_cor = function (ests, axis, R = 2000)
{
    wrap_cor = function (rand)
    {
        beta <- ests$mean$beta
        psiR <- ests$mean$psiR
        psiS <- ests$mean$psiS
        if (rand) {
            ests$se$psiR[is.na(ests$se$psiR)] <- 0
            ests$se$psiS[is.na(ests$se$psiS)] <- 0
            beta <- beta + ests$se$beta * rnorm(length(ests$se$beta), 0, 1)
            psiR <- psiR + ests$se$psiR * rnorm(length(ests$se$psiR), 0, 1)
            psiS <- psiS + ests$se$psiS * rnorm(length(ests$se$psiS), 0, 1)
        }
        expl  <- as.matrix(apply(site[,c("temperature","salinity")], 2, scale))
        coef  <- cbind(beta[,-1], psiR, psiS)
        score <- with(ests$mean, var(cbind(expl, scoreR, scoreS)))
        # Correlation
        coef  <- coef [    ,axis, drop=FALSE]
        score <- score[axis,axis, drop=FALSE]
        cov <- coef %*% score %*% t(coef)
        cor <- cov / sqrt(outer(diag(cov), diag(cov)))
        cor[is.nan(cor)] = 0
        colnames(cor) <- rownames(cor) <- name$Scientific_name
        return(cor)
    }
    cor  <- wrap_cor(FALSE)
    pval <- 0
    for (r in 1:R) pval = pval + (wrap_cor(TRUE) < 0)
    pval <- 2 * pmin(pval, R - pval) / R
    idx  <- order(eigen(cor)$vectors[,1])
    list(mean = cor, pval = pval, idx = idx)
}
cor <- calc_cor(ests, axis=1:7)
sigcor <- cor$mean
sigcor[cor$pval > 0.05] <- 0
sigcor <- sigcor[cor$idx, cor$idx]

png("fig/Fig_E2_corrplot.png")
corrplot::corrplot(
    sigcor, diag=FALSE, type="lower", method="square", tl.cex=0.2, tl.srt=45,
    tl.col=1, addgrid.col=NA)
dev.off()
rm(sigcor, calc_cor)

#------------------------------------------------------------------------------#
# Coefficent plots
coef <- with(ests$mean, cbind(beta[,-1], psiR, psiS))
coef_se <- with(ests$se, 1.96 * cbind(beta[,-1], psiR, psiS))
coef_se[is.na(coef_se)] <- 0
coef[coef_se>abs(coef)] <- NA
coef <- coef[cor$idx,]
species <- name$Scientific_name[cor$idx]

axis <- c("Temperature","Salinity", paste("Local LV",1:ng), paste("Regional LV",1:nf))
df <- lapply(1:ncol(coef), function(k) data.frame(species=species, axis=axis[k], val=coef[,k]))
df <- do.call(rbind, df)
df$species <- factor(df$species, levels=species)
df$axis    <- factor(df$axis   , levels=axis)
df <- df[!is.na(df$val),]

gp_3 <- ggplot(df, aes(x=axis, y=species)) +
    geom_raster(aes(fill=val), interpolate=FALSE) +
    scale_fill_gradient2(low="deeppink", mid="beige", high="deepskyblue", midpoint=0) +
    labs(fill = "Value") +
    theme_bw() + xlab("GLLVM Coefficients") + ylab("Species") +
    theme(
        axis.text.x = element_text(angle=45, hjust=1),
        axis.text.y = element_text(size=4))

ggsave("fig/Fig_S4_coefplot.png", gp_3, width=6, height=12)
rm(cor, df, coef, coef_se, species, axis)

#------------------------------------------------------------------------------#
# Interquiantile of niche conservatism
name$Family <- as.factor(name$Family)
family_set <- lapply(levels(name$Family), function(x) {
    id = which(name$Family == x)
    list(family = x, id = id, sp = name$Scientific_name[id], nsp = length(id))
})
family_set <- family_set[sapply(family_set, function(x) x$nsp) > 9]

# Standardized coefficients
pars <- c("Temperature","Salinity",sprintf("Local LV %s",1:ng),sprintf("Regional LV %s",1:nf))
coef <- data.frame(with(ests$mean, cbind(beta[,-1], psiR, psiS)))
colnames(coef) = pars
var  <- apply(coef^2, 2, sum) / (nS + 1 - c(1,1,1:ng,1:nf))
coef <- coef / matrix(sqrt(var), nrow(coef), ncol(coef), byrow=TRUE)

df <- lapply(family_set, function(set) {
    out = data.frame(
        Family = set$family,
        value  = unlist(apply(coef[set$id,], 2, IQR)),
        Factor = factor(pars, levels=pars)
    )
    rownames(out) <- NULL
    out
})

gp_4 = ggplot(do.call(rbind, df), aes(x=Factor, y=value)) +
    geom_boxplot() +
    geom_point(aes(color=Family, group=Family)) +
    theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x = "", y = "IQR of Standardized Species Response")

gp_5 = ggplot(do.call(rbind, df), aes(x=Family, y=value)) +
    geom_boxplot() +
    geom_point(aes(color=Factor, group=Factor)) +
    theme_bw() + theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x = "", y = "IQR of Standardized Species Response")

ggsave("fig/Fig_3_IQR_lv.png" , gp_4, width=6, height=6)
ggsave("fig/Fig_E3_IQR_fm.png", gp_5, width=6, height=6)

# End