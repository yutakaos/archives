#' Plot posterior densities of SIMMs
#' 
#' \code{plot_post} draws posterior densities from ordinary and beta-dependent SIMMs.
#'
#' @param simm_outs
#'    A simm (and mcmc.list) object, whose posterior is normal.
#' @param simm_outs_beta
#'    A simm (and mcmc.list) object, whose posterior is beta.
#' @param prob
#'    A numeric scalar in the interval (0,1), giving highest posterior density interval. 
#' @param type
#'    The output type. The default is "source", where only source contributions are used
#'    in posterior density plots.
#' 
#' @return
#' A ggplot output.
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
#' 
#' # Posterior density plot
#' plot_post(out_L, out_H, type = "source")

plot_post = function (
    simm_outs, simm_outs_beta = NULL, prob = 0.95, type = c("source", "raw"))
{
    type = match.arg(type)
    num_sources = attr(out_L, "outpar")[1]
    only_L = is.null(simm_outs_beta)
    out_L = as.matrix(simm_outs)
    out_H = if (only_L) out_L else as.matrix(simm_outs_beta)
    
    names = colnames(out_L)
    if (!all(names == colnames(out_H))) stop("Different column names between outputs.")
    if (type == "source") {
        out_L = out_L[,1:num_sources]
        out_H = out_H[,1:num_sources]
    }
    if ("lp" %in% names) {
        idx = which(names == "lp")
        out_L = out_L[,-idx]
        out_H = out_H[,-idx]
        names = names[-idx]
    }
    
    if (only_L) hpd = matrix(0, ncol(out_L), 2)
    else        hpd = coda::HPDinterval(coda::mcmc(out_H), prob)
    df = lapply(1:ncol(out_L), function(k) {
        pd = density(out_L[,k], from = 0, to = 1)
        ud = which(hpd[k,1] < pd$x & pd$x < hpd[k,2])  # which is underdetermined
        out = data.frame(
            source  = rep(names[k], length(pd$x)),
            probs   = pd$x,
            density = pd$y)
        if (length(ud) == 0) {
            out = rbind(cbind(out[1,], post = "beta"), cbind(out, post = "ordinary"))
            out$density[1] = 0
        }
        else {
            out = rbind(cbind(out, post = "ordinary"), cbind(out[ud,], post = "beta"))
            out$density[ud] = 0
        }
        out$post = factor(out$post, levels = c("beta","ordinary"))
        out
    })
    df = do.call(rbind, df)
    df$source = factor(df$source, levels = names)
    
    ggplot(df, aes_string(x = "probs", y = "density", fill = "post")) +
        geom_area() +
        facet_wrap(source ~ ., nrow = 2) + 
        xlim(c(0,1)) +
        guides(fill = "none") +
        xlab("Relative contribution") +
        ylab("Posterior density") + 
        theme_bw() +
        theme(axis.text.y = element_blank(), axis.ticks = element_blank())
}

# End