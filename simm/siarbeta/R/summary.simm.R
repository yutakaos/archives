#' Summary statistics for Bayesian stable isotope mixing models
#' 
#' \code{summary.simm} produces summary statistics (mean, SD, quantile, Gelman-Rbin diagnostic
#' and mean determination coefficient).
#' 
#' @param object
#'    A simm (and mcmc.list) object, which has more than one chain.
#' @param quantiles
#'    A vector of quantiles to evaluate credible intervals.
#' @param ...
#'    Further arguments.
#' 
#' @return
#' A data.frame where each row represents summary statistic of SIMM outputs:
#' \tabular{ll}{
#' \code{Mean}    \tab \code{:} posterior mean \cr
#' \code{SD}      \tab \code{:} posterior standard deviation \cr
#' \code{2.5}     \tab \code{:} posterior 2.5 percent quantile \cr
#' \code{50}      \tab \code{:} posterior median \cr
#' \code{97.5}    \tab \code{:} posterior 97.5 percent quantile \cr
#' \code{rhat}    \tab \code{:} Gelman-Rubin convergence diagnostic \cr
#' \code{mean_r2} \tab \code{:} Mean determination coefficient \cr
#' }
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
#' # Summary statistics
#' summary.simm(out_L)

summary.simm = function (object, quantiles = c(0.025, 0.5, 0.975), ...)
{
    if(!inherits(object, "simm")) stop("object should be 'simm' class.")
    
    mcmc_matrix = coda::mcmc(as.matrix(object))
    outs = summary(mcmc_matrix, quantiles, ...)  # coda:::summary.mcmc
    outs = cbind(outs$statistics[,c("Mean", "SD")], outs$quantiles)
    
    rhat = coda::gelman.diag(object, autoburnin = FALSE, multivariate = FALSE)
    outs = cbind(outs, rhat = pmax(1, rhat$psrf[,1]))
    
    outpar = attr(object, "outpar")
    num_sources = outpar[1]
    mean_r2 = rep(NA, sum(outpar))
    r2 = cor(mcmc_matrix[,1:num_sources])^2
    mean_r2[1:num_sources] = (apply(r2, 1, sum) - 1) / (num_sources - 1)
    outs = cbind(outs, mean_r2 = mean_r2)
    
    outs
}

# End