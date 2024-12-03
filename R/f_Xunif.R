#' Standardized effects for a Uniform covariate
#'
#' `f_Xunif()` returns a list to build various standardized effects (specified through the `model` argument) under the assumption that the corresponding covaiate follows a Uniform distribution (either, discrete or continuous). The first three elements of the list are designed to be used in INLA to specify an effect in its standardized version through the model `generic0`. The list also reports the scaling constant, a sample from the Uniform distribution on the covariate, and  the basis matrix evaluated at this sample: these elements can be useful for the user. The argument `plot_check=TRUE` returns a plot showing realizations of the standardized effect along with the histograms of the mean and variance of each realization.

#' @param x A numeric vector (for `model="linear","pspline1","pspline2","pspline_2D"`) or an ordered factor with more than 1 level (for `model="iid","rw1","rw2","besag"``).
#' @param model String to select the type of effect.
#' @param K A positive integer larger than 3, specifying the number of basis functions. Needed for `model="pspline1","pspline2","pspline_2D"`.
#' @param adj_mat Matrix needed for `model="besag"`, representing the adjacency matrix from the lattice.
#' @param m A number indicating the lower boundary of the support of X. If not provided, it is set to the minimum value from `x`. Needed for `model="linear","pspline1","pspline2","pspline_2D"`.
#' @param M A number indicating the upper boundary of the support of X. If not provided, it is set to the maximum value from `x`.  Needed for `model="linear","pspline1","pspline2","pspline_2D"`.
#' @param fixed Logical, indicating whether the effect must be treated as fixed (TRUE) or random (FALSE). If it is treated as fixed, a 0 mean constraint is imposed By default, TRUE. The only model affected is `model="iid"`.
#' @param plot_check Logical, indicating whether a plot to check the standardization must be printed (TRUE) or not (FALSE). By default, FALSE.
#' @param n_sim Positive integer, indicating the number of realizations to generate for the check plot.
#'
#' @examples
#' K <- 20
#' x <- factor(sample(1:K, 100, replace = TRUE), ordered = TRUE, levels = c(1:K))
#' f_Xunif(x, model = "iid", plot_check = TRUE)
#' f_Xunif(x, model = "iid", fixed = FALSE, plot_check = TRUE)
#' f_Xunif(x, model = "rw1", plot_check = TRUE)
#' f_Xunif(x, model = "rw2", plot_check = TRUE)
#'
#' K <- 20
#' x <- runif(1000, min = 2, max = 5)
#' f_Xunif(x, model = "linear", plot_check = TRUE)
#' f_Xunif(x, model = "pspline1", K = K, plot_check = TRUE)
#' f_Xunif(x, model = "pspline2", K = K, plot_check = TRUE)
f_Xunif <- function(
    x, model, K,
    adj_mat,
    m = NULL, M = NULL, scale_Q=T,
    fixed = TRUE, plot_check = F, n_sim = 1000) {
  # Discrete ------
  if (model == "iid") {
    result <- iid_standard(x = x, fixed = fixed)
  }
  if (model == "rw1") {
    result <- rw_standard(x = x, order = 1)
  }
  if (model == "rw2") {
    result <- rw_standard(x = x, order = 2)
  }
  if (model == "besag") {
    result <- besag_standard(x = x, adj_mat = adj_mat)
  }
  # Continuous -----
  if (model == "linear") {
    result <- linear_standard(x = x, m = m, M = M)
  }
  if (model == "pspline1") {
    result <- pspline_standard(x = x, K = K, order = 1, m = m, M = M)
  }
  if (model == "pspline2") {
    result <- pspline_standard(x = x, K = K, order = 2, m = m, M = M)
  }
  if (model == "pspline_2D") {
    result <- Pspline_2D_standard(x = x, K = K, m = m, M = M)
  }

  # Precision matrox unscaled -----
  if (!scale_Q) {
    result$precision <- result$precision/result$scaling_constant
  }

  # Check results through simulation-----

  if (plot_check) {

    if (model!="pspline_2D" & model !="besag") {
      realizations <- result$basis_distribution %*%
        t(mvtnorm::rmvnorm(n_sim, sigma = gen_inv(result$precision)))

    return(
      ggpubr::annotate_figure(
        ggpubr::ggarrange(
          ggplot2::ggplot(
            data =
              tidyr::gather(as.data.frame(realizations))
          ) +
            ggplot2::geom_line(alpha = 0.1, ggplot2::aes(
              x = rep(result$X_distribution, n_sim),
              y = value, group = key
            )) +
            ggplot2::theme_light() +
            ggplot2::labs(y = "f(x)", x = "x"),
          check_GMRF(realizations, fixed = fixed,X_distribution = result$X_distribution),
          ncol = 1
        ),
        top = paste("Model:", model)
      )
    )
    } else {
      return(
        ggpubr::annotate_figure(
          check_GMRF(realizations, fixed = fixed,X_distribution = result$X_distribution),
          top = paste("Model:", model)
        )
      )
    }
  }

  return(result)

}
