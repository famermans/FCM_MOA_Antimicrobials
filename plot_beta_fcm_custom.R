plot_beta_fcm_custom <- function (x, color = NA, shape = NA, size = NA, labels = c("Factor 1",
                                                                                   "Factor 2",
                                                                                   "Factor 3"), legend.pres = NULL) 
{
  legend.ops <- NULL
  if (sum(is.na(color)) > 0) 
    color = rep("f1", nrow(x$points))
  if (sum(is.na(shape)) > 0) {
    shape = rep("f2", nrow(x$points))
    legend.ops <- FALSE
  }
  if (sum(is.na(size)) > 0) {
    size = rep("f3", nrow(x$points))
    legend.ops <- FALSE
  }
  var.pcoa <- vegan::eigenvals(x)/sum(vegan::eigenvals(x))
  PcoA <- as.data.frame(x$points)
  names(PcoA)[1:2] <- c("Axis1", "Axis2")
  PcoA <- cbind(PcoA, color, shape, size)
  ggplot2::ggplot(PcoA, ggplot2::aes(x = Axis1, y = Axis2,
                                     color = color, shape = shape, size = size)) + ggplot2::geom_point(alpha = 0.7,
                                                                                                       size = 4) +
    ggplot2::geom_point(colour = "grey90", size = 1.5) + 
    ggplot2::scale_color_manual(values = c("#a65628", "red", 
                                           "#ffae19", "#4daf4a", "#1919ff", "darkorchid3", 
                                           "magenta")) +
    ggplot2::labs(x = paste0("Axis1 (", 
                             round(100 * var.pcoa[1], 1), "%)"), y = paste0("Axis2 (", 
                                                                            round(100 * var.pcoa[2], 1), "%)")) +
    ggplot2::ggtitle("Ordination of phenotypic fingerprints") +
    ggplot2::labs(color = labels[1], shape = labels[2], size = size[3]) +
    ggplot2::guides(color = legend.pres, shape = legend.ops, size = legend.pres)
}
