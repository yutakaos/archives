#------------------------------------------------------------------------------#
# Multispecies nonlinear Ricker model
#------------------------------------------------------------------------------#

make_data_multisp_ricker = function (
    tl = 100, nsp = 40, connectance = 0.2, burnin = 200)
{
    generate_M = function (nsp, connectance, x) {
        count <- 0
        while(count < 1000) {
            M <- -rbinom(nsp*nsp, 1, connectance)*runif(nsp*nsp, 0.05, 0.1)
            M <- matrix(M, nrow=nsp)
            diag(M) <- 0
            diag(M) <- -(1 + M %*% x)/x
            # check IGR > 0
            igr <- sapply(1:nsp, function(i) {
                z <- rep(0,nsp)
                z[-i] <- solve(M[-i,-i]) %*% rep(-1,nsp-1)
                1 + M[i,] %*% z
            })
            if (all(igr > 0)) break
            count <- count + 1
        }
        if (count==1000) stop("Could not find interaction matrix.")
        return(M)
    }
        
    # make data
    r <- 3  # growth parameter
    while(1) {
        x <- runif(nsp, 0.5, 1)  # biomass
        M <- generate_M(nsp, connectance, x)
        out <- matrix(NA, burnin+tl, nsp)
        out[1,] <- x + rnorm(nsp, 0, 1e-2)
        for (t in 2:(burnin+tl)) {
            out[t,] <- out[t-1,] * exp(r*(1 + M%*%out[t-1,]))
            out[t, out[t,]<1e-6] <- 0
            out[t, out[t,]>1e+1] <- 10
        }
        coexist_sp <- out[burnin+tl,]!=0 & out[burnin+tl,]!=10
        #if(nsp*0.9 <= sum(coexist_sp)) break
        if(nsp == sum(coexist_sp)) break
    }
    M <- M[coexist_sp, coexist_sp]
    out <- out[1:tl+burnin, coexist_sp]
    colnames(out) <- paste0("sp.",1:ncol(out))
    list(df=data.frame(t = 1:tl, out), M=M, nsp=sum(coexist_sp))
}

#------------------------------------------------------------------------------#
# Multispecies VAR model
#------------------------------------------------------------------------------#

make_data_multisp_var = function (
    tl = 100, nsp = 40, connectance = 0.2, burnin = 200)
{
    generate_M = function (nsp, connectance) {
        count <- 0
        while(count < 1000) {
            M <- rbinom(nsp*nsp,1,connectance)*runif(nsp*nsp, 0.1, 0.3)
            M <- M*(2*rbinom(nsp*nsp,1,0.5)-1)  # sign
            M <- matrix(M, nrow=nsp)
            diag(M) <- 0.1
            if (abs(eigen(M)$values)[1]<1) break
            count <- count + 1
        }
        if (count==1000) stop("Could not find interaction matrix.")
        return(M)
    }
    M <- generate_M(nsp, connectance)
    out <- matrix(NA, burnin+tl, nsp)
    out[1,] <-  rnorm(nsp, 0, 1e-2)
    for (t in 2:(burnin+tl)) {
        out[t,] <- M %*% out[t-1,] + 0.5 * rnorm(nsp)
    }
    out <- out[1:tl+burnin,]
    colnames(out) <- paste0("sp.",1:nsp)
    list(df=data.frame(t = 1:tl, out), M=M, nsp=nsp)
}

# End