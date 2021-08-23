#' Underdetermination diagnositc for Bayesian stable isotope mixing models
#' 
#' \code{evaluate_ump} produces underdetermination diagnostics.
#' 
#' @param simm_outs
#'    A simm (and mcmc.list) object, whose posterior is beta.
#' @param prob
#'    A numeric scalar in the interval (0,1), giving highest posterior density (HPD) interval. 
#' 
#' @return
#' A data.frame where each column represents:
#' \tabular{ll}{
#' \code{lower} \tab \code{:} lower bound of HPD interval of beta-dependent SIMM \cr
#' \code{upper} \tab \code{:} upper bound of HPD interval of beta-dependent SIMM \cr
#' \code{delta} \tab \code{:} the difference between lower and upper bounds \cr
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
#' # Run MCMC for beta-dependent SIMM
#' out_H = siarbeta(
#'     mixture, sources, correct = 0, concdep = 1, alpha = 1, beta = 1000,
#'     error = "variance", chains = 3, iters = 10000, burns = 5000, thins = 4)
#' 
#' evaluate_ump(out_H, prob = 0.95)

evaluate_ump = function (simm_outs, prob = 0.95)
{
    if (attr(simm_outs, "beta") < 200) {
        beta = attr(simm_outs, "beta")
        if (beta > 1) {
            message = paste0(
                "Beta may be too small to evaluate underdetermination (beta = ",
                beta, "). Our recommendation is beta >= 200.")
            warning(message)
        }
        else {
            message = paste0("simm_outs is ordinary SIMM outputs (beta = ", beta, ").")
            stop(message)
        }
    }
    mcmc_matrix = coda::mcmc(as.matrix(simm_outs))
    hpdi  = coda::HPDinterval(mcmc_matrix, prob)
    delta = hpdi[,2] - hpdi[,1]
    delta[length(delta)] = NA
    outs = cbind(hpdi, delta = delta)
    outs
}

# End