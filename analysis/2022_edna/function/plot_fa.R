#--------------------------------------------------------------------------------#
# Plot for Factor Analysis
#--------------------------------------------------------------------------------#

.plot_fa = function (mat, main = "", n = 10) {
    k  <- 1 + floor(n * (mat / max(abs(range(mat))) + 1))
    pl <- cm.colors(1 + 2 * n)
    plot.new()
    plot.window(xlim=c(0,ncol(mat)), ylim=c(0,nrow(mat)), xaxt="n", yaxt="n")
    for (i in 1:ncol(mat)) for(j in 1:nrow(mat)) {
        polygon(c(i-1,i,i,i-1), c(j-1,j-1,j,j), col=pl[k[j,i]], border=NA)
    }
    polygon(c(0,ncol(mat),ncol(mat),0), c(0,0,nrow(mat),nrow(mat)))
    title(main)
}

.plot_fa_varimax = function (psi, score, rott) {
    nf <- ncol(psi)
    if (nf == 0) {
        plot.new()
        plot.new()
    }
    else if (nf == 1 || !rott) {
        .plot_fa(score, main = "scores (location)")
        .plot_fa(psi  , main = "loadings (species)")
    }
    else {
        varimax = varimax(psi)
        psiv   = varimax$loadings[,1:nf]
        scorev = score %*% varimax$rotmat
        .plot_fa(scorev, main = "scores (location)")
        .plot_fa(psiv  , main = "loadings (species)")
    }
}

plot_report = function (report, record, rott = FALSE)
{
    if (!is.null(dev.list()["RStudioGD"])) dev.off(dev.list()["RStudioGD"])
    par(mfrow = c(3,2))
    psiQ   <- with(report, cbind(1, psiR))
    scoreQ <- with(report, cbind(alpha, scoreR))
    with(report, .plot_fa_varimax(psiS, scoreS, rott=rott))
    with(report, .plot_fa_varimax(psiQ, scoreQ, rott=rott))
    
    if (length(report$kappa) != 0) {
        barplot(sqrt(8)/report$kappa, col="grey", main="Range (rho)")
    }
    else plot.new()
    
    if ((nr <- nrow(record)) != 1)
    {
        idx  <- max(2, nr-20):nr
        iter <- record$iter[idx]
        nll_ <- record$nll [idx]
        nll_[is.infinite(nll_)] <- NA
        plot(iter, nll_, type="l", ylab="nll", main="Negative log-liklihood")
    }
}

# End
