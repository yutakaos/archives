#--------------------------------------------------------------------------------#
# Make inputs for TMB model
#--------------------------------------------------------------------------------#

make_tmb_inputs = function (data, nF, nG)
{
    ip <- NULL
    nS <- ncol(data$cm)  # number of species
    nL <- nrow(data$cm)  # number of locations
    nE <- ncol(data$expl)
    nM <- data$spde$n.spde
    
    ## data
    ip$data <- list(
        Y = data$cm,
        X = data$expl,
        Penalty = c(0, 0, 0, 0),
        G0 = data$spde$param.inla$M0,
        G1 = data$spde$param.inla$M1,
        G2 = data$spde$param.inla$M2
    )
    
    ## parameters
    ip$pars <- list(
        beta  = matrix(0, nS, nE),
        psiSv = rep(0.1, nS * nF - nF * (nF - 1) / 2),
        psiRv = rep(0.1, nS * nG - nG * (nG - 1) / 2),
        alpha = rep(0, nL),
        omegaS = matrix(0, nM, nF),
        scoreR = matrix(0, nL, nG),
        log_sigma = 0,
        log_kappa = rep(0, nF)
    )
    
    # Fix parameters
    ip$map <- list()
    
    # Declare random
    ip$random <- c("alpha", "omegaS", "scoreR")
    
    ## return
    return(ip)
}

# End
