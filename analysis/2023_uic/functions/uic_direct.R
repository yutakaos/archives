#------------------------------------------------------------------------------#
# uic_direct
#------------------------------------------------------------------------------#

uic_direct = function (
    block, lib = c(1, NROW(block)), pred = lib, group = NULL,
    lib_var = 1, tar_var = 2, cond_var = NULL,
    norm = 2, E = 1, tau = 1, tp = 0, nn = "e+1", num_surr = 1000, alpha = 0.05, 
    exclusion_radius = NULL, epsilon = NULL,
    is_naive = FALSE, knn_method = c("KD","BF"))
{
    X <- tar_var
    X <- X[!X %in% c(lib_var,cond_var)]
    
    # Wrapper function
    uic_wrap = function (block, tar_var, cond_var = NULL, tp = 0)
    {
        uic.optimal(
            block, lib, pred, group, lib_var, tar_var, cond_var,
            norm, E, tau, tp, nn, num_surr, alpha,
            exclusion_radius, epsilon, is_naive, knn_method)
    }
    
    # Make delay time-series
    make_dc = function (X, lib = c(1,NROW(X)), E = 0, group=NULL)
    {
        lib <- rbind(lib)
        if (is.null(group)) {
            X <- as.matrix(X)
            NAsF <- matrix(NA, nrow=pmax(0, E), ncol=ncol(X))
            NAsB <- matrix(NA, nrow=pmax(0,-E), ncol=ncol(X))
            nNA <- nrow(NAsB)
            out <- lapply(1:nrow(lib), function(i) {
                L1 <- lib[i,1]
                L2 <- lib[i,2]
                dc <- rbind(NAsF, X[L1:L2,,drop=FALSE], NAsB)
                dc [nNA+1:(L2-L1+1),,drop=FALSE]
            })
            out <- do.call(rbind, out)
            return(out)
        }
        lib <- lapply(1:nrow(lib), function(i) {
            L <- lapply(unique(group), function(g) range((lib[i,1]:lib[i,2])[group==g]))
            L <- do.call(rbind, L)
        })
        lib <- do.call(rbind, lib)
        make_dc(X, lib, E)
    }
    
    # 1st Round: X -> Y(lib_var)
    out <- lapply(X, function(x) 
        data.frame(uic_wrap(block, x, cond_var, tp), target=x, direct="NA")
    )
    out <- do.call(rbind, out)
    out <- out[order(-out$te),]
    out$direct <- factor(out$direct, levels=c("+","-","NA"))
    out$direct[out$pval < alpha] <- "+"
    
    # 2nd Round: X -> Z -> Y(lib_var)
    Lag = NULL
    for (i in 1:(nrow(out)-1)) {
        if (out$direct[i] != "+") next
        lag <- make_dc(block[out$target[i]], lib, -out$tp[i], group)
        colnames(lag) <- sprintf(".lag.%s", i)
        Lag <- cbind(Lag, lag)
        for (j in (i+1):nrow(out))
        {
            if (out$direct[j] != "+") next
            Cond <- c(cond_var, colnames(Lag))
            outj <- with(out, 
                uic_wrap(data.frame(block, Lag), target[j], Cond, tp[j]))
            if (outj$pval < alpha) next
            out$direct[j] <- "-"
        }
    }
    out[order(as.numeric(out$direct)),]
}

# End