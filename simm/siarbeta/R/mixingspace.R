#' Plot isotopic mixing space
#' 
#' \code{mixingspace} draws isotopic mixing space.
#' 
#' @param mixture
#'    A matrix with each target individual as a seperate row and each isotope as
#'    a seperate column.
#' @param sources
#'    A matrix containing the mean and standard deviations of the source values for
#'    each of the isotopes.
#' @param axis
#'    The element axes used as x and y in mixing space.
#' @param source_names
#'    The source names.
#' @param element_names
#'    The element names.
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
#' mixingspace(
#'     mixture, sources,
#'     source_names  = paste("Source", 1:nrow(sources)),
#'     element_names = c("d13C", "d15N"))

mixingspace = function (
    mixture, sources, axis = 1:2, source_names = NULL, element_names = NULL)
{
    axis = axis[1:2]
    mixture = as.data.frame(mixture)
    sources = as.data.frame(sources)
    
    mixture = mixture[,axis]
    sources = sources[,c(sapply(axis, function(k) 2*(k-1)+1:2))]
    
    if (is.null(source_names )) source_names  = paste("Source" , 1:nrow(sources))
    if (is.null(element_names)) element_names = paste("Element", 1:2)
    
    colnames(mixture) = c("E1", "E2")
    colnames(sources) = c("E1", "E1e", "E2", "E2e")
    sources = cbind(sources, xmin = sources$E1 - sources$E1e)
    sources = cbind(sources, xmax = sources$E1 + sources$E1e)
    sources = cbind(sources, ymin = sources$E2 - sources$E2e)
    sources = cbind(sources, ymax = sources$E2 + sources$E2e)
    sources = cbind(sources, Sources = source_names)
    
    ggplot(sources, aes_string(x = "E1", y = "E2", color = "Sources")) +
        geom_point() +
        geom_point(aes_string(x = "E1", y = "E2"), data = mixture, color = "black") +
        geom_errorbarh(aes_string(xmin = "xmin", xmax = "xmax"), height = 0) +
        geom_errorbar (aes_string(ymin = "ymin", ymax = "ymax"), width  = 0) +
        xlab(element_names[1]) +
        ylab(element_names[2]) +
        theme_bw() +
        theme(legend.position = "bottom", legend.direction = "vertical") +
        theme(legend.background = element_rect(color = "black"))
}

# End