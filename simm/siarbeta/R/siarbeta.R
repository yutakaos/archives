#' Beta-dependent stable isotope mixing model
#' 
#' Applications of beta-dependent stable isotope mixing model.
#' \code{siarbeta} diagnoses the uncertainty due to underdetermination by beta-posteriors.
#'
#' @param mixture
#'    A matrix with each target individual as a seperate row and each isotope as a seperate column.
#' @param sources
#'    A matrix containing the mean and standard deviations of the source values for
#'    each of the isotopes.
#' @param correct
#'    A matrix containing the mean and standard deviations of the fractional correction
#'    values for each of the isotopes.
#' @param concdep
#'    A matrix containing the mean and standard deviations of the concentration dependence
#'    values for each of the isotopes. Note that this function does not use the standard deviations.
#' @param alpha
#'    A scalar or vector containing the Dirichlet prior parameters.
#' @param beta
#'    The beta in beta-dependent posterior probability. The default is 1 (ordinary SIAR).
#' @param error
#'    The error structure. The default is "stock".
#' @param chains
#'    The number of chains. The default is 3.
#' @param iters
#'    The number of iterations to run. The default is 5000.
#' @param burns
#'    The size of the burn-in. The default is 1000.
#' @param thins
#'    The amount of thinning of the iterations. The default is 4.
#' @param silent
#'    prevents messages from being printed to the R console. The default is False.
#' @param source_names
#'    The source names.
#' 
#' @return
#' A mcmc.list output.
#' 
#' @examples
#' mixture = matrix(0, nrow = 1, ncol = 2)
#' sources = matrix(c(
#'     -1, 0.25,  1, 0.25,
#'     -1, 0.25, -1, 0.25,
#'      1, 0.25,  1, 0.25,
#'      1, 0.25, -1, 0.25), nrow = 4, byrow = TRUE)
#' 
#' # Run MCMC for ordinary SIMM
#' out_L = siarbeta(
#'     mixture, sources, correct = 0, concdep = 1, alpha = 1, beta = 1,
#'     error = "variance", chains = 3, iters = 10000, burns = 5000, thins = 4)
#' 
#' # Run MCMC for beta-dependent SIMM
#' out_H = siarbeta(
#'     mixture, sources, correct = 0, concdep = 1, alpha = 1, beta = 1000,
#'     error = "variance", chains = 3, iters = 10000, burns = 5000, thins = 4)

siarbeta = function (
    mixture, sources, correct = 0, concdep = 1, alpha = 1, beta = 1,
    error = c("stock", "variance", "residual", "parnell"),
    chains = 3, iters = 5000, burns = 1000, thins = 4, silent = FALSE,
    source_names = NULL)
{
    error = match.arg(error)
    num_sources  = nrow(sources)
    num_elements = ncol(mixture)
    
    mixture = as.matrix(mixture)
    sources = as.matrix(sources)
    correct = as.matrix(correct)
    concdep = as.matrix(concdep)
    
    if (length(correct) == 1) correct = matrix(correct, num_sources, 2*num_elements)
    if (length(concdep) == 1) concdep = matrix(concdep, num_sources, 2*num_elements)
    if (length(alpha)   == 1) alpha   = rep(alpha, num_sources)
    if (is.null(source_names)) source_names = paste("Source", 1:num_sources)
    
    if (!silent) {
        cat("\nBayesian Stable Isotope Mixing Models ")
        cat(paste0("(SIAR, beta = ", beta, ")\n"))
    }
    
    error = switch(error, variance = 0, residual = 1, parnell = 2, stock = 3)
    if (nrow(mixture) == 1 && error != 0) {
        if (!silent) warning("'variance' error is used for single mixture sample.\n")
        error = 0
    }
    
    model = new(SIMM)
    model$set_srcX_summary(sources)
    model$set_srcD_summary(correct)
    model$set_srcQ_summary(concdep)
    model$set_mixX(mixture)
    model$set_priors(alpha)
    model$set_temp(beta, beta)
    
    valid_output = "ps"
    if (error == 1) valid_output = c(valid_output, "sd")
    if (error == 2) valid_output = c(valid_output, "sd")
    if (error == 3) valid_output = c(valid_output, "ep")
    valid_output = c(valid_output, "lp")
    
    t0 = proc.time()[3]
    howmany = floor(iters / 20)
    outs = rep(list(NULL), chains)
    for (ch in 1:chains) {
        model$init(error)
        out = list(ps = NULL, sd = NULL, ep = NULL, lp = NULL)
        if (!silent) cat(paste("Chain", ch, ": |"))
        for (i in 1:iters) {
            model$update(1)
            if (burns < i && (i - burns) %% thins == 0) {
                x = model$output(TRUE)
                out$ps = rbind(out$ps, x$ps)
                out$sd = rbind(out$sd, x$sd)
                out$ep = rbind(out$ep, x$ep)
                out$lp = rbind(out$lp, x$lp)
            }
            if (!silent && i %% howmany == 0) cat(ifelse(i <= burns, "==", "++"))
        }
        if (!silent) cat("|\n")
        colnames(out$ps) = source_names
        colnames(out$sd) = paste0("SD"     , 1:num_elements)
        colnames(out$ep) = paste0("epsilon", 1:num_elements)
        colnames(out$lp) = "lp"
        out = coda::mcmc(do.call(cbind, out[valid_output]))
        attr(out, "mcpar") = c(burns + thins, iters, thins)
        outs[[ch]] = out
    }
    outs = coda::mcmc.list(outs)
    
    if (!silent) {
        cat("Job completed successfully.\n")
        cat(paste("Elapsed time:", round(proc.time()[3] - t0, 2), "seconds.\n"))
    }
    
    outpar = c(num_sources, 0, 0, 1)
    if (error == 1) outpar[2] = num_elements  # residual
    if (error == 2) outpar[2] = num_elements  # parnell
    if (error == 3) outpar[3] = num_elements  # stock
    
    attr(outs, "outpar") = outpar
    attr(outs, "beta"  ) = beta
    attr(outs, "class" ) = c("simm","mcmc.list")
    outs
}

# End