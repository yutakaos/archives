#------------------------------------------------------------------------------#
# Run Generalized Linear Latent Variable Model
#------------------------------------------------------------------------------#
# Set directory
setwd("D:/Desktop/2022_edna")
if (!dir.exists("fig")) dir.create("fig")
if (!dir.exists("out")) dir.create("out")

# Load library
library(geosphere)
library(INLA)
library(ggplot2)
library(cowplot)
library(maps)
theme_set(theme_cowplot())

# Read data
name <- read.csv("data/20211110_species.csv", header=TRUE)
site <- read.csv("data/20211110_site.csv"   , header=TRUE)
cm <- read.csv("data/20211110_community_matrix.csv", header=TRUE)
dc <- read.csv("data/20211110_costdist.csv", header=TRUE)
dc <- as.matrix(dc)

# Map Polygon
Jpn <- map_data("world", region="Japan")

#------------------------------------------------------------------------------#
# Model-input data
#  * Remove sites with missing environmental data
#  * Remove species complex taxa
#  * Remove data-poor species (Occurrence <= 3 sites)
idxS <- name$input == 1
idxL <- site$input == 1
name <- name[idxS,]
site <- site[idxL,]
cm <- cm[idxS, idxL]
dc <- dc[idxL, idxL]

# Scale explanatory data
expl <- cbind(1, apply(site[,c("temperature","salinity")], 2, scale))

# Geographic distance
dg <- geosphere::distm(site[,c('lon','lat')]) / 1000

# Make INLA mesh
remove_inner_mesh = function (mesh, dg, dc, cutoff = 1, err = 30)
{
    nl = nrow(dg)
    tv = mesh$graph$tv
    # Indices of inner meshes
    op = NULL
    for (i in 1:(nl-1)) for (j in (i+1):nl) {
        if (dg[i,j] == 0) next
        if (dc[i,j] - err < cutoff * dg[i,j]) next
        vs  = mesh$idx$loc[c(i,j)]
        idx = which(sapply(apply(tv, 1, intersect, y = vs), length) == 2)
        op  = c(op, idx)
    }
    inn = sort(unique(op))
    # Remove inner meshes
    mesh$graph$tv  = mesh$graph$tv [-inn,]
    mesh$graph$tt  = mesh$graph$tt [-inn,]
    mesh$graph$tti = mesh$graph$tti[-inn,]
    mesh
}
mesh <- inla.mesh.create(site[,c('lon','lat')])
mesh <- remove_inner_mesh(mesh, dg, dc)
#plot(mesh, draw.vertices = TRUE)

# SPDE with Matern distance
spde <- inla.spde2.matern(mesh, alpha=2)

# Adjacent Map
loc <- with(mesh, data.frame(lon=loc[,1], lat=loc[,2]))
tv <- with(mesh$graph, rbind(tv[,1:2], tv[,2:3], tv[,c(3,1)]))
tv <- t(apply(tv, 1, sort))
tv <- tv[!duplicated(1e5*tv[,1]+tv[,2]),]
tv <- data.frame(
    xi = loc[tv[,1],1], xe = loc[tv[,2],1], 
    yi = loc[tv[,1],2], ye = loc[tv[,2],2]
)

gp_1 = ggplot(loc, aes(x=lon, y=lat)) +
    coord_fixed(xlim=range(loc$lon), ylim=range(loc$lat)) +
    geom_segment(aes(x=xi, y=yi, xend=xe, yend=ye), data=tv) +
    geom_polygon(aes(x=long, y=lat, group=group, fill=region), data=Jpn, fill="grey70") +
    geom_point() +
    theme_bw() + xlab("Longitude") + ylab("Latitude")

ggsave("fig/Fig_S2_AdjacentMap.png", gp_1, width=6, height=6)

# Output R Data
data <- list(
    name = name,
    site = site,
    cm = 1 * t(as.matrix(cm) != 0),
    expl = expl,
    mesh = mesh,
    spde = spde
)
rm(list=ls()[ls() != "data"])
#save(data, file="out/edna_japan_data.RData")

#------------------------------------------------------------------------------#
#load('out/edna_japan_data.RData')

# Load library
library(TMB)

# Load functions
source("function/make_tmb_inputs_lvm.R")
source("function/plot_fa.R")

# Run GLLVM
nf <- 3      # number of latent variables (local)
ng <- 2      # number of latent variables (regional)

# Link DLL
compile("function/src/edna_lvm.cpp")
dyn.load(dynlib("function/src/edna_lvm"))

# Initialization
filename <- sprintf("op_lvm_f%s_g%s.RData", nf, ng)
ip <- make_tmb_inputs(data, nF=nf, nG=ng)
ip$data$Penalty <- c(100, 0, 0, 0)

# Make TMB object
obj <- with(ip, MakeADFun(data=data, parameters=pars, random=random, map=map, hessian=FALSE))
obj$env$inner.control <- list(step.tol=1e-6, tol10=1e-6, grad.tol=1e-6, maxit=100)

# first run
init <- NULL
init$fn <- obj$fn(obj$par)
init$gr <- obj$gr(obj$par)
save.image(paste0("out/", filename))

# run model
iter <- 0
record <- data.frame(iter=iter, nll=NA)
while (1) {
    opt <- nlminb(
        start = obj$env$last.par[-obj$env$random],
        objective = obj$fn, gradient = obj$gr,
        control = list(eval.max=50, iter.max=1e3, trace=0, rel.tol=1e-6))
    iter <- iter + opt$iterations
    record <- rbind(record, c(iter=iter, nll=opt$objective))
    save.image(paste0("out/", filename))
    plot_report(obj$report(obj$env$last.par), record, rott=FALSE)
    if (opt$convergence == 0 || opt$message == "false convergence (8)")
        break
}
print(opt$message); print(opt$objective)
# Try optim instead of nlminb if opt$message == "false convergence (8)".

# Hessian matrix
est  <- obj$report(obj$env$last.par.best)
hess <- optimHess(
    obj$env$last.par.best[-c(obj$env$random)], obj$fn, obj$gr,
    control=list(maxit=100, reltol=1e-6))

# Summary
summary_outs = function (opt_outs)
{
    # Covariance
    est  <- opt_outs$est
    hess <- opt_outs$hess
    cov  <- MASS::ginv(as.matrix(hess))
    se_all <- sqrt(diag(abs(cov)))
    
    # Parameters
    pars <- colnames(hess)
    nl <- nrow(est$scoreS)
    ns <- nrow(est$beta)
    ne <- ncol(est$beta)
    nf <- ncol(est$scoreS)
    ng <- ncol(est$scoreR)
    
    # SE
    psiSt <- matrix(NA, nf, ns)
    psiRt <- matrix(NA, ng, ns)
    psiSt[upper.tri(psiSt, diag=TRUE)] <- se_all[pars=="psiSv"]
    psiRt[upper.tri(psiRt, diag=TRUE)] <- se_all[pars=="psiRv"]
    se <- list(
        beta = matrix(se_all[pars=="beta"], ns, ne),
        log_kappa = se_all[pars=="log_kappa"],
        log_sigma = se_all[pars=="log_sigma"],
        psiS = t(psiSt),
        psiR = t(psiRt)
    )
    
    # Confidence intervals
    lower <- list(
        beta  = est$beta - 1.96 * se$beta,
        kappa = exp(log(est$kappa) - 1.96 * se$log_kappa),
        sigma = exp(log(est$sigma) - 1.96 * se$log_sigma),
        psiS  = est$psiS - 1.96 * se$psiS,
        psiR  = est$psiR - 1.96 * se$psiR
    )
    upper <- list(
        beta  = est$beta + 1.96 * se$beta,
        kappa = exp(log(est$kappa) + 1.96 * se$log_kappa),
        sigma = exp(log(est$sigma) + 1.96 * se$log_sigma),
        psiS  = est$psiS + 1.96 * se$psiS,
        psiR  = est$psiR + 1.96 * se$psiR
    )
    
    # Outputs
    list(mean = est, se = se, lower = lower, upper = upper)
}
opt_outs <- list(obj=obj, est=est, hess=hess)
ests <- summary_outs(opt_outs)

save(opt_outs, ests, filename, nf, ng, 
     file = sprintf("op_lvm_f%s_g%s_ests.RData", nf, ng))

# End