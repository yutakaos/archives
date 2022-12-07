#------------------------------------------------------------------------------#
# GLLVM analysis 2
#------------------------------------------------------------------------------#
# Set directory
setwd("D:/Desktop/2022_edna")
if (!dir.exists("fig")) dir.create("fig")

# Load library
library(rfishbase)
library(polycor)
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())

# Load data
name <- read.csv("data/20211110_species.csv", header=TRUE)
name <- name[name$input == 1,]

# Load RData
load("out/op_lvm_f3_g2_ests.RData")
rm(nf, ng, opt_outs, filename)

#------------------------------------------------------------------------------#
# Ecological information from FishBase
fb <- data.frame(rfishbase::ecology(name$Scientific_name))
# Whether species include ecological data
idx <- which(!is.na(fb$SpecCode) & !is.na(fb$EcologyRefNo))
# Whether ecological data are registered for at least 40 fish species
fb <- fb[,apply(fb==-1, 2, sum, na.rm=TRUE) >= 40]

# Polychoric correlation
pcor <- sapply(names(fb), function(i) sapply(names(fb), function(j) polychor(fb[idx,i], fb[idx,j])))
rownames(pcor) <- colnames(pcor) = colnames(fb)
print(pcor)
# Remove ecological data which is correlated to other data
fb <- fb[,!colnames(fb) %in% c("Mangroves","Sand","Mud","Rocky")]

#--------------------------------------------------------------------------------#
ggplot_glm = function (df)
{
    df$var = factor(df$var, levels=names(fb))
    ggplot(df, aes(x = x, y = y)) +
        geom_point (color="darkblue") +
        geom_jitter(height=0.01) +
        stat_smooth(data = df[ df$sig,], method="glm", formula = y ~ x,
                    method.args=list(family="binomial"), color="darkblue") +
        stat_smooth(data = df[!df$sig,], method="glm", formula = y ~ x,
                    method.args=list(family="binomial"), color="darkblue", linetype="dotted") +
        scale_y_continuous(breaks=c(0,1)) +
        facet_wrap(. ~ var, nrow=3) +
        labs(x="Species-specific response", y="Habitat type") +
        theme_bw()
}

# GLM
glm_out <- list(
    R1 = summary(glm(y ~ ., data = data.frame(y = ests$mean$psiR[idx,1], -fb[idx,]))),
    R2 = summary(glm(y ~ ., data = data.frame(y = ests$mean$psiR[idx,2], -fb[idx,])))
)
df1 <- lapply(names(fb), function(x)  # for local LV 1
    data.frame(
        x = ests$mean$psiR[idx,1], y = -fb[idx,x], var = x,
        sig = glm_out$R1$coefficients[x,4] < 0.05)
)
df2 <- lapply(names(fb), function(x)  # for local LV 2
    data.frame(
        x = ests$mean$psiR[idx,2], y = -fb[idx,x], var = x,
        sig = glm_out$R2$coefficients[x,4] < 0.05)
)
df1 <- do.call(rbind, df1)
df2 <- do.call(rbind, df2)
gp_1 <- ggplot_glm(df1)
gp_2 <- ggplot_glm(df2)

ggsave("fig/Fig_S5_fb_llv1.png", gp_1, width=6, height=6)
ggsave("fig/Fig_S6_fb_llv2.png", gp_2, width=6, height=6)

# End