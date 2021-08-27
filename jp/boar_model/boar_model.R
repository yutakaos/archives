#--------------------------------------------------------------------------------#
# R script for boar model
#--------------------------------------------------------------------------------#
# Load environmental data
data <- read.csv("data/envi.csv", header = TRUE)

# Habitat area
area <- data[,"all_forest"] + data[,"abandoned"]

# Environmental variables to explain growth rates
env_var_names <- c(
   "brl_forest",   # broad-leaf forests
   "forest_edge",  # forest edges
   "bamboo",       # bamboo forests
   "field",        # agricultural fields
   "orchard",      # orchards
   "abandoned",    # abandoned fields
   "urban")        # urban lands
env_vars <- asin(data[,env_var_names] / data[,"area"])
env_vars <- apply(env_vars, 2, scale)  # standardization

#--------------------------------------------------------------------------------#
# Load trap data
data <- read.csv("data/trap.csv", header = TRUE)
trap08 <- data[data["year"] == 2008,]
trap09 <- data[data["year"] == 2009,]

# Box trap
capt <- cbind(trap08[,"captured"], trap09[,"captured"])
trap <- cbind(trap08[,"trapspan"], trap09[,"trapspan"])

idcs <- apply(!is.na(capt * trap), 2, which)
ncs  <- sapply(idcs, length)
idcs <- sapply(idcs, function(x) c(x, rep(NA, max(ncs) - length(x))))

# Number of hunted individuals
hunt <- cbind(trap08[,"hunted"], trap09[,"hunted"])
capt[is.na(capt)] <- 0
hunt <- hunt + capt

#--------------------------------------------------------------------------------#
# Weights for spatial autocorrelations
weis <- read.csv("data/border_length.csv",header = TRUE)
weis <- weis[,-(1:2)]
weis[is.na(weis)] <- 0
weis <- t(sapply(weis, function(x) c(x) / sum(x + 1e-6)))
diag(weis) <- -1

idss <- which(apply(weis, 1, sum) != -1)
nss  <- length(idss)

#--------------------------------------------------------------------------------#
# JAGS inputs
ns <- nrow(env_vars)
nv <- ncol(env_vars)

inputs <- list(
   # input data
   data <- list (
      NS   = ns,    # number of cities
      NV   = nv,    # number of variabels
      NCS  = ncs,   # number of cities with observations
      IDCS = idcs,  # city IDs with observations
      CAPT = capt,  # number of captured boars (by box traps)
      TRAP = trap,  # number of box traps
      HUNT = hunt,  # number of hunted boars
      AREA = area,  # habitat area
      ENVI = cbind(icpt = 1, env_vars),  # environmental variables
      NSS  = nss,   # number of non-isolated cities
      IDSS = idss,  # non-isolated city IDs 
      WEIS = weis   # weights of auto-correlations
   ),
   
   # intializations (vs and g_coef are undeclared)
   inits <- function() list(
      r_capt = 0.5,
      num_0  = hunt[,1] + 10,
      num_g  = hunt + 10,
      sigma  = 1,
      coef   = rep(0, nv),
      g_rand = rep(0, ns),
      g_prec = 0.1,
      s_prec = 0.1
   ),
   
   # parameters to save
   params <- c(
      "r_capt","density","growth","g_coef","g_rand",
      "sigma" ,"g_prec" ,"s_prec")
)

#--------------------------------------------------------------------------------#
library(R2jags)
post.bugs <- with(inputs, jags(
   data, inits, params, "boar_model.txt",
   n.chains=3, n.iter=2000000, n.burnin=1000000, n.thin=1000))

# Outputs
post.bugs$BUGSoutput

# g_coef
g_coef <- data.frame(post.bugs$BUGSoutput$sims.matrix[,paste0("g_coef[",1:8,"]")])
colnames(g_coef) <- c("intercept", env_var_names)

outs <- lapply(g_coef, function(b) {
   vs <- b != 0
   c(quantile(b[vs], probs = c(0.5, 0.025, 0.975)), pr = mean(vs))
})
outs <- do.call(rbind, outs)
outs

# End