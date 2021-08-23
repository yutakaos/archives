#' Correlation plot
#' 
#' \code{plot_corr} draws probability density plots (diagonal), heat maps (lower) and
#' correlation values (upper).
#'
#' @param simm_outs
#'    A simm (and mcmc.list) object.
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
#' # Correlation plot
#' plot_corr(out_L[,1:nrow(sources)])

plot_corr = function (simm_outs)
{
    color_palette = colorRampPalette(c("white", "lightblue", "yellow", "red"))
    lower_list = list(
        continuous = function (data, mapping) 
        {
            ggplot(data = data, mapping = mapping) + 
                stat_density2d(aes_string(fill="..density.."), geom="tile", contour = FALSE) +
                scale_fill_gradientn(colors = color_palette(100))
        }
    )
    upper_list = list(
        continuous = wrap(ggally_cor, stars = FALSE)
    )
    
    outs = as.data.frame(as.matrix(simm_outs))
    ggpairs(outs, lower = lower_list, upper = upper_list) +
        theme_test()
}

# End