#------------------------------------------------------------------------------#
# CCM test function
#------------------------------------------------------------------------------#

ccm_test = function (df, lib_var, tar_var, maxE=10, nsample=100)
{
    E <- rEDM::EmbedDimension(
        dataFrame=df, lib=c(1,nrow(df)), pred=c(1,nrow(df)),
        maxE=maxE, Tp=1, columns=lib_var, target=lib_var,
        showPlot=FALSE, noTime=TRUE, numThreads=1)
    E <- E$E[which.max(E$rho)]
    
    out_ccm <- lapply(setdiff(tar_var, lib_var), function(x) {    
        # surrogate test
        ts0 <- rEDM::SurrogateData(df[[x]], method="ebisuzaki", num_surr=nsample)
        out_surr <- lapply(1:nsample, function(i) {
            with(rEDM::CCM(
                dataFrame=cbind(df, surr=ts0[,i]),
                E=E, Tp=-1, columns=lib_var, target="surr",
                libSizes=nrow(df)-E+1, random=FALSE, includeData=TRUE),
            CCM1_PredictStat)
        })
        out_surr <- do.call(rbind, out_surr)
        # CCM
        out <- with(rEDM::CCM(
            dataFrame=df, E=E, Tp=-1, columns=lib_var, target=x,
            libSizes=nrow(df)-E+1, random=FALSE, includeData=TRUE),
            CCM1_PredictStat)
        out_min <- with(rEDM::CCM(
            dataFrame=df, E=E, Tp=-1, columns=lib_var, target=x,
            libSizes=E+5, sample=nsample, includeData=TRUE),
            CCM1_PredictStat)
        pval1 <- mean(out_min $RMSE < out$RMSE)
        pval2 <- mean(out_surr$RMSE < out$RMSE)
        # output
        data.frame(
            lib=lib_var, tar=x, out, RMSE0=mean(out_surr$RMSE),
            pval=max(pval1,pval2))
    })
    do.call(rbind, out_ccm)
}

#------------------------------------------------------------------------------#
# TE test function
#------------------------------------------------------------------------------#

te_test = function (df, lib_var, tar_var, nsample=100)
{
    out_te <- lapply(setdiff(tar_var, lib_var), function(x) {    
        out <- RTransferEntropy::transfer_entropy(
            df[,x], df[,lib_var], lx=1, ly=1, quiet=TRUE,
            shuffles=1, nboot=nsample)
        data.frame(lib=lib_var, tar=x, te=out$coef[1,1], pval=out$coef[1,4])
    })
    do.call(rbind, out_te)
}

#------------------------------------------------------------------------------#
# Utility function for multispeices test
#------------------------------------------------------------------------------#

make_M_from_test <- function (out_test, nsp, pval="pval") {
    M0 <- matrix(1, nrow=nsp, ncol=nsp)
    diag(M0) <- 0
    M0[M0==1] <- out_test[,pval]
    t(M0)
}

calc_metrics = function (M, M0) {
    # from pcalg package (version 2.7-12)
    compute_shd = function(m1, m2) {
        shd <- 0
        s1 <- m1 + t(m1)
        s2 <- m2 + t(m2)
        s1[s1 == 2] <- 1
        s2[s2 == 2] <- 1
        ds <- s1 - s2
        ind <- which(ds > 0)
        m1[ind] <- 0
        shd <- shd + length(ind)/2
        ind <- which(ds < 0)
        m1[ind] <- m2[ind]
        shd <- shd + length(ind)/2
        d <- abs(m1 - m2)
        shd + sum((d + t(d)) > 0)/2
    }
    
    # metrics
    nsp <- nrow(M)
    df_test <- data.frame(data=1*(M[diag(nsp)==0]!=0), pred=M0[diag(nsp)==0])
    data.frame(
        tpr = with(df_test, sum(pred[data==1]< 0.05)/sum(data==1)),  # recall
        tnr = with(df_test, sum(pred[data==0]>=0.05)/sum(data==0)),  # speficity
        acc = with(df_test, sum((data==1)-(pred<0.05)==0)/nrow(df_test)),  # accuracy
        shd = compute_shd(1*(M!=0), 1*(M0<0.05))  # structural Hamming distance
    )
}

# End