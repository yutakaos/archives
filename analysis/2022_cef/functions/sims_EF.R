#------------------------------------------------------------------------------#
# Simulate biomass at equilibrium
#------------------------------------------------------------------------------#

sims_EF = function (R = NA, S = 30, c = 0.2, p = 0.5, apr = FALSE)
{
    #------------------------------------------------------------#
    # R   : dummy variables
    # S   : number of species
    # c   : connectance
    # p   : proportion of positive interactions
    # apr : whether approximate biomasses are computed
    #------------------------------------------------------------#
    # Expectations and variances of demographic parameters
    Er = 1.5; Vr = 1/3  # intrinsic growth rate
    Es = 6.0; Vs = 1/3  # self-regulation
    
    # Intervals for uniform distributions 
    Ir <- sqrt(12*Vr)/2
    Is <- sqrt(12*Vs)/2
    
    # Non-diagonal element ID
    id <- matrix(1:(S*S), S, S)
    id <- id[lower.tri(id) | upper.tri(id)]
    
    # Simulate biomass at equilibrium (x) 
    LA <- S * (S - 1) * c *  p 
    LB <- S * (S - 1) * c * (1 - p)
    LA <- floor(LA) + rbinom(1, 1, LA - floor(LA))
    LB <- floor(LB) + rbinom(1, 1, LB - floor(LB))
    while (1)
    {
        si <- runif(S, Es-Is, Es+Is)
        A <- B <- matrix(0, S, S)
        A[sample(id, LA)] <- runif(LA, 0, 1)
        B[sample(id, LB)] <- runif(LB, 0, 1)
        M  <- diag(si, S) - A + B
        iM <- try(solve(M), silent = TRUE)
        if (!"try-error" %in% class(iM)) break
    }
    ri <- runif(S, Er-Ir, Er+Ir)
    xi <- c(iM %*% ri)
    if (!apr) return(list(x = xi))
    
    ## Compute approximate biomass (y)
    uAi <- apply(A, 1, sum)
    uBi <- apply(B, 1, sum)
    uAo <- apply(A, 2, sum)
    uBo <- apply(B, 2, sum)
    sumA <- sum(uAo)
    sumB <- sum(uBo)
    
    rR  <- sum(uAo * ri)  / sumA; rC  <- sum(uBo * ri)  / sumB
    sR  <- sum(uAo * si)  / sumA; sC  <- sum(uBo * si)  / sumB
    bRR <- sum(uAo * uAi) / sumA; bCR <- sum(uBo * uAi) / sumB
    bRC <- sum(uAo * uBi) / sumA; bCC <- sum(uBo * uBi) / sumB
    
    H  <- (sR - bRR) * (sC + bCC) + bCR*bRC
    R  <- (rR*sC + rR*bCC - rC*bRC) / H
    C  <- (rC*sR + rR*bCR - rC*bRR) / H
    yi <- (ri + uAi*R - uBi*C) / si
    list(x = xi, y = yi)
}

# End