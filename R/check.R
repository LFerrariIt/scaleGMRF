#' Plots to check standardization of a GMRF.
#'
#' `check_GMRF()` returns a plot made up by 1 or more panels to check the standardization of a Gaussian effect. If a sample of the distribution of the covariate is specified, a plot of the realizations is also returned.
#'
#' @param realizations A numeric matrix, containing a realization of the process in each column.
#' @param fixed Logical, indicating whether the effect must be treated as fixed or random. By default, `TRUE`.
#' @param X_dist A numeric matrix, containing values sampled from the distribution of X. By default, `NULL`.
#'
#' @return Plot.
#'
#' @details This function is used to visually check the standardization of the Gaussian effect. It is called in the functions `standardize_GMRF()` and `standardize_X_unif()` of the `scaleGMRF` package when `plot_check=T`. The plot provided in the output is different for fixed and random effects: in both cases, the "Mean" value in the plot should be close to 1 if standardization has been correctly applied. See `vignette("standardization",package="scaleGMRF")` for more details.
#'
check_GMRF <- function(realizations, fixed = T, X_dist = NULL) {
  if (!is.matrix(realizations)) {
    stop("`realizations` must be a matrix object. Instead, it is a ", class(realizations), " object.")
  }

  if (ncol(realizations) < 2 | nrow(realizations) < 0) {
    stop("`realizations` must have at least 2 columns and 2 rows. Instead, it has dimension ", nrow(realizations), "x", ncol(realizations), ".")
  }

  if (any(is.na(realizations))) {
    stop("`realizations` must not contain NA values.")
  }

  if (!is.logical(fixed)) {
    stop("fixed must be logical.")
  }

  n_sim <- ncol(realizations)

  if (fixed) {
    # mean_sample <- colMeans(realizations)

    var_sample <- apply(realizations, 2, function(x) var(x) * (length(x) - 1) / length(x))

    plot_1 <-
      # ggpubr::ggarrange(
      # ggplot2::ggplot() +
      # ggplot2::geom_histogram(ggplot2::aes(mean_sample), fill = "grey") +
      # ggplot2::geom_vline(xintercept = 0, ggplot2::aes(col = "0")) +
      # ggplot2::geom_vline(
      # xintercept = mean(mean_sample), lty = "dashed",
      # ggplot2::aes(col = "Mean")) +
      # ggplot2::theme_light() +
      # ggplot2::labs(title = base::expression(E[X] ~ "[" ~ f(X) ~ "]"), x = paste("Mean:", round(mean(mean_sample), 2))),
      ggplot2::ggplot() +
      ggplot2::geom_histogram(ggplot2::aes(var_sample), fill = "grey") +
      ggplot2::geom_vline(
        xintercept = mean(var_sample), lty = "dashed"
      ) +
      ggplot2::theme_light() +
      ggplot2::labs(title = base::expression(Var["i=1"]^n ~ "[" ~ f(x[i]) ~ "]"), y = "Number of realizations", x = paste("Within realization variance\n Mean:", round(mean(var_sample), 2)))
    # )
  } else {
    var_sample <- apply(realizations, 1, var)
    if (is.null(X_dist)) {
      plot_1 <- ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = 1:length(var_sample), y = var_sample)) +
        ggplot2::geom_hline(yintercept = mean(var_sample), lty = "dashed") +
        ggplot2::theme_light() +
        ggplot2::labs(title = base::expression(Var["b=1"]^B ~ "[" ~ f^"(b)" ~ (x) ~ "]"), x = paste("x\n Mean:", round(mean(var_sample), 2)), y = "Between realizations variance")
    } else {
      plot_1 <- ggplot2::ggplot() +
        ggplot2::geom_line(ggplot2::aes(x = as.matrix(X_dist)[, 1], y = var_sample)) +
        ggplot2::geom_hline(yintercept = mean(var_sample), lty = "dashed") +
        ggplot2::theme_light() +
        ggplot2::labs(title = base::expression(Var["b=1"]^B ~ "[" ~ f^"(b)" ~ (x) ~ "]"), x = paste("x\n Mean:", round(mean(var_sample), 2)), y = "Between realizations variance")
    }
  }

  if (is.null(X_dist)) {
    return(plot_1 +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)))
  } else {
    return(
      ggpubr::ggarrange(
        ggplot2::ggplot(
          data =
            tidyr::gather(as.data.frame(realizations))
        ) +
          ggplot2::geom_line(alpha = 0.1, ggplot2::aes(
            x = rep(as.matrix(X_dist)[, 1], n_sim),
            y = value, group = key
          )) +
          ggplot2::theme_light() +
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
          ggplot2::labs(y = "f(x)", x = "x"),
        plot_1 + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)),
        ncol = 2
      )
    )
  }
}
