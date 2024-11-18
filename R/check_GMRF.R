#' Checks if the effect has been properly standardized.
#'
#' `check_GMRF()` returns a plot with three panels. Realizations of the effect from the argument are generated and displayed in the first panel. The expected value and variances of the realizations with respect to the Uniform distribution of the covariate are computed and plotted via histograms. The means of both statistics are also reported, which should be equal respectively to 0 (only if `fixed=TRUE`) and 1.
#'
#' @param result A list returned from the `f_Xunif()` function.
#'
#' @return Plot.
#'
check_GMRF <- function(result) {
  n_sim <- 1000
  f_sample <- result$basis_distribution %*%
    t(mvtnorm::rmvnorm(n_sim,
      sigma = gen_inv(result$precision)
    ))

  mean_sample <- colMeans(f_sample)

  var_sample <- apply(f_sample, 2, var)

  return(
    ggpubr::ggarrange(
      ggplot2::ggplot(
        data =
          tidyr::gather(as.data.frame(f_sample))
      ) +
        ggplot2::geom_line(alpha = 0.1, ggplot2::aes(
          x = rep(as.matrix(result$X_distribution)[, 1], n_sim),
          y = value, group = key
        )) +
        ggplot2::theme_light() +
        ggplot2::labs(y = "f(x)", x = "x"),
      ggpubr::ggarrange(
        ggplot2::ggplot() +
          ggplot2::geom_histogram(ggplot2::aes(mean_sample), fill = "grey") +
          ggplot2::geom_vline(xintercept = 0) +
          ggplot2::geom_vline(xintercept = mean(mean_sample), col = 2) +
          ggplot2::theme_light() +
          ggplot2::labs(x = "Mean", title = round(mean(mean_sample), 2)),
        ggplot2::ggplot() +
          ggplot2::geom_histogram(ggplot2::aes(var_sample), fill = "grey") +
          ggplot2::geom_vline(xintercept = 1) +
          ggplot2::geom_vline(xintercept = mean(var_sample), col = 2) +
          ggplot2::theme_light() +
          ggplot2::labs(x = "Variance", title = round(mean(var_sample), 2))
      ),
      ncol = 1
    )
  )
}
