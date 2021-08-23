#--------------------------------------------------------------------------------#
# Demo
#--------------------------------------------------------------------------------#

library(siarbeta)
library(ggplot2)
library(gridExtra)
dir.create("figures")

#--------------------------------------------------------------------------------#
# Example 1: Toy simulations
#--------------------------------------------------------------------------------#

mixture = matrix(0, nrow = 1, ncol = 2)
sources = matrix(c(
    -1, 0.25,  1, 0.25,
    -1, 0.25, -1, 0.25,
     1, 0.25,  1, 0.25,
     1, 0.25, -1, 0.25), nrow = 4, byrow = TRUE)

# Run MCMC (L = ordinary, H = beta)
set.seed(0204)
out_L = siarbeta(
    mixture, sources, correct = 0, concdep = 1, alpha = 1, beta = 1,
    error = "variance", chains = 3, iters = 10000, burns = 5000, thins = 4)

out_H = siarbeta(
    mixture, sources, correct = 0, concdep = 1, alpha = 1, beta = 1000,
    error = "variance", chains = 3, iters = 10000, burns = 5000, thins = 4)

# Summary statistics
summary(out_L)
summary(out_H)
evaluate_ump(out_H)

# Plot mixing space
g1 = mixingspace(
    mixture, sources,
    source_names  = paste("Source", 1:4),
    element_names = c("d13C", "d15N") )

# Correlation plot
g2 = plot_corr(out_L[,1:4])

# Posterior density plot
g3 = plot_post(out_L, out_H, type = "source")

# Output figures
ggsave("figures/Fig_1.png", plot = grid.arrange(g1, g3, ncol = 2), width = 8, height = 6)
ggsave("figures/Fig_2.png", plot = g2, width = 8, height = 8)

#--------------------------------------------------------------------------------#
# Example 2: Geese field data
#--------------------------------------------------------------------------------#

library(siar)
data("geese1demo", "sourcesdemo", "correctionsdemo", "concdepdemo")
mixture = geese1demo
sources = sourcesdemo[,-1]
correct = correctionsdemo[,-1]
concdep = concdepdemo[,-1]
source_names = as.character(sourcesdemo[,1])
rm(geese1demo, sourcesdemo, correctionsdemo, concdepdemo)

# Run MCMC (L = ordinary, H = beta)
set.seed(0204)
out_L = siarbeta(
    mixture, sources, correct, concdep, alpha = 1, beta = 1,
    error = "parnell", source_names = source_names,
    chains = 3, iters = 5000, burns = 1000, thins = 4)

out_H = siarbeta(
    mixture, sources, correct, concdep, alpha = 1, beta = 1000,
    error = "parnell", source_names = source_names,
    chains = 3, iters = 5000, burns = 1000, thins = 4)

# Summary statistics
summary(out_L)
summary(out_H)
evaluate_ump(out_H)

# Plot mixing space
g1 = mixingspace(
    mixture, sources,
    source_names  = source_names,
    element_names = c("d13C", "d15N") )

# Correlation plot
g2 = plot_corr(out_L[,1:4])

# Posterior density plot
g3 = plot_post(out_L, out_H, type = "source")

# Output figures
ggsave("figures/Fig_3.png", plot = grid.arrange(g1, g3, ncol = 2), width = 8, height = 6)
ggsave("figures/Fig_4.png", plot = g2, width = 8, height = 8)

# End