#' Checks if the effect has been properly standardized.
#'
#' `check_GMRF()` returns a plot with 2 panels: the means and variances of the realizations are plotted via histograms. The means of both statistics are also reported, which should be equal respectively to 0 (only if `fixed=TRUE`) and 1.
#'
#' @param realizations A numeric matrix, containing a realization of the process in each column.
#'
#' @return Plot.
#'
check_GMRF <- function(realizations, fixed = T,X_distribution=NULL) {
  if (!is.matrix(realizations)) {
    stop("`realizations` must be a matrix.")
  }

  if (ncol(realizations) < 2 | nrow(realizations) < 0) {
    stop("`realizations` must have at least 2 columns and 2 rows. Instead, it has dimension ", nrow(x), "x", ncol(x), ".")
  }

  if (any(is.na(realizations))) {
    stop("`realizations` must not contain NA values.")
  }

  if (fixed) {
    mean_sample <- colMeans(realizations)

    var_sample <- apply(realizations, 2, function(x) var(x) * (length(x) - 1) / length(x))

    return(
      ggpubr::ggarrange(
        ggplot2::ggplot() +
          ggplot2::geom_histogram(ggplot2::aes(mean_sample), fill = "grey") +
          ggplot2::geom_vline(xintercept = 0, ggplot2::aes(col = "0")) +
          ggplot2::geom_vline(
            xintercept = mean(mean_sample), lty = "dashed",
            ggplot2::aes(col = "Mean")
          ) +
          ggplot2::theme_light() +
          ggplot2::labs(title = base::expression(E[X]~"["~f(X)~"]"), x = paste("Mean:", round(mean(mean_sample), 2))),
        ggplot2::ggplot() +
          ggplot2::geom_histogram(ggplot2::aes(var_sample), fill = "grey") +
          ggplot2::geom_vline(
            xintercept = mean(var_sample), lty = "dashed") +
          ggplot2::theme_light() +
          ggplot2::labs(title = base::expression(Var[X]~"["~f(X)~"]"), x = paste("Mean:", round(mean(var_sample), 2)))
      )
    )
  } else {
    if (is.null(X_distribution)) {
      var_sample <- apply(realizations, 1, var)

      return(
        ggplot2::ggplot() +
          ggplot2::geom_line(ggplot2::aes(x = 1:length(var_sample), y = var_sample)) +
          ggplot2::geom_hline(yintercept = mean(var_sample), lty = "dashed") +
          ggplot2::theme_light() +
          ggplot2::labs(title="Variance",y=base::expression(Var["i=1"]^n~"["~f[i](x)~"]"), x = "x")
      )
    } else {
      var_sample <- apply(realizations, 1, var)

      return(
        ggplot2::ggplot() +
          ggplot2::geom_line(ggplot2::aes(x = as.matrix(X_distribution)[,1], y = var_sample)) +
          ggplot2::geom_hline(yintercept = mean(var_sample), lty = "dashed") +
          ggplot2::theme_light() +
          ggplot2::labs(title="Variance",y=base::expression(Var["i=1"]^n~"["~f[i](x)~"]"), x = "x")
      )
    }
  }
}
