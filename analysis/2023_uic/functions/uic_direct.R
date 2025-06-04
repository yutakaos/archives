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
    
    # 1st Round: X -> Y(lib_var)
    out <- lapply(X, function(tar_var) {
        out.x <- uic.optimal(
            block, lib, pred, group, lib_var, tar_var, cond_var,
            norm, E, tau, tp, nn, num_surr, alpha, TRUE,
            exclusion_radius, epsilon, is_naive, knn_method)
        data.frame(out.x, target=tar_var)
    })
    out <- do.call(rbind, out)
    seq_test <- rep(0, nrow(out))
    idx <- which(out$seq_test>0)
    seq_test[idx[which.max(out$ete[idx])]] <- 1
    
    # 2nd Round: X -> Z -> Y(lib_var)
    lag <- block[0]
    while (1) {
        idx <- which(out$seq_test>0 & seq_test==0)
        if (length(idx) == 0) break
        z <- which.max(seq_test)
        lag <- cbind(lag, rUIC::make_block(block, out$target[z], -out$tp[z], group))
        out.cond <- lapply(idx, function(j)
            out.cond <- uic.optimal(
                cbind(block,lag), lib, pred, group,
                lib_var, tar_var=out$target[j], c(cond_var,colnames(lag)),
                norm, E, tau, tp=out$tp[j], nn, num_surr, alpha, FALSE,
                exclusion_radius, epsilon, is_naive, knn_method)
        )
        out.cond <- do.call(rbind, out.cond)
        r <- max(seq_test)
        seq_test[idx[which.max(out.cond$ete)]] <- r+1
        seq_test[idx[out.cond$pval>=alpha]]    <- -r
    }
    out$seq_test <- seq_test
    return(out)
}


# compute stability
fun_smap = function (out)
{
    #out <- out_1[out_1$pval<0.05 ,]
    theta <- c(
        0, 1e-04, 3e-04, 0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 0.5,
        0.75, 1, 1.5, 2, 3, 4, 6, 8)
    
    smap <- lapply(1:nrow(taxa), function(i) {
        # make block
        outi <- out[out$effect==i,]
        lib_var <- taxa$Taxa_ID[i]
        B <- matrix(NA, nrow=nrow(df), ncol=0)
        cond_var <- NULL
        if (nrow(outi)!=0) {
            B <- lapply(1:nrow(outi), function(k) {
                x <- taxa$Taxa_ID[outi$cause[k]]
                rUIC::make_block(df, x, -outi$tp[k]-1, "plot")
            })
            B <- data.frame(do.call(cbind, B))
            cond_var <- colnames(B)
        }
        
        # simplex
        simp <- rUIC::simplex(
            data.frame(plot=df$plot, effect=df[,lib_var], B), group="plot",
            lib_var="effect", cond_var=cond_var, E=1:10, tp=1, alpha=0.05)
        optE <- with(simp, max(0, E[pval < 0.05]))
        
        # smap
        lib <- lapply(unique(df$plot), function(x) range(which(df$plot==x)))
        lib <- do.call(rbind, lib)
        L <- lapply(0:optE-1, function(k) rUIC::make_block(df, lib_var, k, "plot"))
        L <- data.frame(do.call(cbind, L))
        DF <- cbind(L,B)
        smap <- lapply(theta, function(x)
            rEDM::SMap(dataFram=DF, lib=lib, pred=lib, E=1*(ncol(DF)==2), Tp=1, theta=x, 
                 columns=colnames(DF)[-1], target=colnames(DF)[1],
                 embedded=(ncol(DF)!=2), noTime=TRUE)
        )
        fit <- sapply(smap, function(X) with(X$predictions, rEDM::ComputeError(Observations, Predictions)))
        smap <- smap[[which.min(unlist(fit["RMSE",]))]]
        coef <- smap$coefficients[-1]
        outi <- cbind(outi, smap=apply(coef[,-1:-ncol(L),drop=FALSE], 2, median, na.rm=TRUE))
        smry <- data.frame(
            cause  = c(rep(lib_var, optE+1), taxa$Taxa_ID[outi$cause]),
            effect = lib_var, E = c(0:optE, -outi$tp))
        list(stats=outi, summary=smry, coef=coef)
    })
    
    # dynamical stability
    smry <- do.call(rbind, lapply(smap, function(x) x$summary))
    coef <- do.call(cbind, lapply(smap, function(x) x$coef))
    coef <- coef[,smry$E!=0]
    smry <- smry[smry$E!=0,]
    smry$i <- match(smry$cause , taxa$Taxa_ID)
    smry$j <- match(smry$effect, taxa$Taxa_ID)
    time <- which(complete.cases(coef))
    stab <- data.frame(df[,1:2], val=NA)
    stab[time,]$val <- sapply(time, function(t) {
        nsp <- nrow(taxa)
        adj <- data.frame(smry, val=unlist(coef[t,]))
        A <- matrix(0, nrow=nsp, ncol=max(adj$E)*nsp)
        for (k in 1:nrow(adj)) {
            A[adj$i[k],(adj$E[k]-1)*nsp+adj$j[k]] <- adj$val[k]
        }
        A <- rbind(A, diag(nrow=ncol(A)-nsp, ncol=ncol(A)))
        abs(eigen(A)$values[1])
    })
    stab = stab[complete.cases(stab),]
    
    # output
    list(stats=do.call(rbind, lapply(smap, function(x) x$stats)), stability=stab)
}


gg_mynet = function (out, name, size=1)
{
    if (length(size)==1) size <- rep(size, length(name))
    mygraph <- igraph::graph_from_data_frame(
        data.frame(from="origin", to=name),
        vertices=data.frame(name=c("origin", name), value=c(NA, size)))
    from <- out$cause  + 1  # +1 (origin)
    to   <- out$effect + 1  # +1 (origin)
    ggraph(mygraph, layout='dendrogram', circular=TRUE) + 
        geom_node_point(aes(filter=leaf, size=value), color="coral") +
        geom_conn_bundle(
            data=get_con(from=from, to=to, col=factor(out$smap>0)),
            aes(x=x*0.95, y=y*0.95, color=col), alpha=0.5, tension=0.7,
            arrow=grid::arrow(angle=20, length=unit(0.1,"inches"), type="closed")) +
        theme_void() + theme(legend.position="none")
}

# End