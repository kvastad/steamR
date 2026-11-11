#' Generate trait score trajectories
#'
#' Plot trait-associated window scores along an ordered cell lineage and
#' calculate a pattern score for each window based on linear and LOESS fits.
#'
#' @param se A Seurat object containing the window-score variables
#'   and the cell-ranking variable in its metadata.
#' @param meta_prefix Character string giving the prefix of the window-score
#'   variables.
#' @param n_vars Positive integer specifying the number of window-score
#'   variables to evaluate.
#' @param rank_var Character string specifying the variable used to order
#'   cells along the lineage. Defaults to "MetaCell_rank".
#' @param ncol Number of columns in the faceted plot. Defaults to 9.
#' @param point_alpha Transparency of the plotted points. Defaults to 0.3.
#' @param point_size Size of the plotted points. Defaults to 0.8.
#' @param smooth_color Color of the LOESS smoothing line. Defaults to "red".
#' @param show_se Logical; whether to display the standard error band.
#'   Defaults to TRUE.
#'
#' @return A ggplot object with one panel per window-score variable.
#'   The pattern scores are stored in the plot's "pattern_scores" attribute.
#'
#' @details
#' For each window-score variable, the function calculates the R-squared
#' of a linear model and a LOESS model. The pattern score is the maximum
#' of the two R-squared values, allowing either a linear or nonlinear
#' model to capture the observed pattern.
#'
#' The pattern score is intended as a relative measure of how well the
#' observed window scores follow a linear or nonlinear pattern along the
#' ordered lineage. It is not a statistical significance test and does
#' not provide a p-value or establish the biological significance of
#' the observed pattern. Higher scores indicate a better fit to the
#' observed data among the two models considered.
#'
#' Negative R-squared values are retained when both models have negative
#' R-squared values. In such cases, the larger (less negative) value is
#' selected as the pattern score.
#'
#' @export
generate_trait_score_trajectory <- function(
    se,
    meta_prefix,
    n_vars,
    rank_var = "MetaCell_rank",
    ncol = 9,
    point_alpha = 0.3,
    point_size = 0.8,
    smooth_color = "red",
    show_se = TRUE
) {
  
  # ----- Input validation -----
  
  if (!is.character(meta_prefix) ||
      length(meta_prefix) != 1L ||
      is.na(meta_prefix) ||
      !nzchar(meta_prefix)) {
    stop("`meta_prefix` must be a non-empty character string.")
  }
  
  if (!is.numeric(n_vars) ||
      length(n_vars) != 1L ||
      is.na(n_vars) ||
      !is.finite(n_vars) ||
      n_vars < 1 ||
      n_vars != floor(n_vars)) {
    stop("`n_vars` must be a positive integer.")
  }
  
  if (!is.character(rank_var) ||
      length(rank_var) != 1L ||
      is.na(rank_var) ||
      !nzchar(rank_var)) {
    stop("`rank_var` must be a non-empty character string.")
  }
  
  if (!is.numeric(ncol) ||
      length(ncol) != 1L ||
      is.na(ncol) ||
      !is.finite(ncol) ||
      ncol < 1 ||
      ncol != floor(ncol)) {
    stop("`ncol` must be a positive integer.")
  }
  
  if (!is.logical(show_se) ||
      length(show_se) != 1L ||
      is.na(show_se)) {
    stop("`show_se` must be TRUE or FALSE.")
  }
  
  # ----- Construct variable names -----
  
  vars <- paste0(meta_prefix, seq_len(n_vars), "_1")
  
  # Fetch the ranking variable and window scores.
  df <- Seurat::FetchData(
    se,
    vars = c(rank_var, vars)
  )
  
  # Confirm that the requested variables are available.
  missing_vars <- setdiff(c(rank_var, vars), colnames(df))
  
  if (length(missing_vars) > 0L) {
    stop(
      "The following variables were not found: ",
      paste(missing_vars, collapse = ", ")
    )
  }
  
  if (!is.numeric(df[[rank_var]])) {
    stop("`rank_var` must contain numeric values.")
  }
  
  # ----- Reshape data -----
  
  df_long <- tidyr::pivot_longer(
    df,
    cols = tidyselect::all_of(vars),
    names_to = "variable",
    values_to = "value"
  )
  
  df_long$variable <- factor(
    df_long$variable,
    levels = vars
  )
  
  # ----- Calculate pattern scores -----
  
  calc_pattern_score <- function(df_sub) {
    
    # Remove missing and non-finite values for model fitting.
    fit_df <- data.frame(
      rank = df_sub[[rank_var]],
      value = df_sub$value
    )
    
    fit_df <- fit_df[
      is.finite(fit_df$rank) & is.finite(fit_df$value),
      ,
      drop = FALSE
    ]
    
    variable_name <- as.character(df_sub$variable[1L])
    
    # Insufficient data or no variation in the response.
    if (nrow(fit_df) < 3L ||
        length(unique(fit_df$rank)) < 2L ||
        length(unique(fit_df$value)) < 2L) {
      return(tibble::tibble(
        variable = variable_name,
        r2_linear = NA_real_,
        r2_loess = NA_real_,
        pattern_score = NA_real_
      ))
    }
    
    # Linear model.
    lm_fit <- tryCatch(
      stats::lm(value ~ rank, data = fit_df),
      error = function(e) NULL
    )
    
    r2_linear <- if (!is.null(lm_fit)) {
      summary(lm_fit)$r.squared
    } else {
      NA_real_
    }
    
    # LOESS model.
    loess_fit <- tryCatch(
      suppressWarnings(
        stats::loess(value ~ rank, data = fit_df)
      ),
      error = function(e) NULL
    )
    
    r2_loess <- NA_real_
    
    if (!is.null(loess_fit)) {
      
      pred_loess <- tryCatch(
        stats::predict(loess_fit, newdata = fit_df),
        error = function(e) NULL
      )
      
      if (!is.null(pred_loess) &&
          length(pred_loess) == nrow(fit_df) &&
          all(is.finite(pred_loess))) {
        
        ss_total <- sum(
          (fit_df$value - mean(fit_df$value))^2
        )
        
        ss_residual <- sum(
          (fit_df$value - pred_loess)^2
        )
        
        if (ss_total > 0) {
          r2_loess <- 1 - ss_residual / ss_total
        }
      }
    }
    
    # Select the higher R-squared value.
    valid_scores <- c(r2_linear, r2_loess)
    valid_scores <- valid_scores[is.finite(valid_scores)]
    
    pattern_score <- if (length(valid_scores) > 0L) {
      max(valid_scores)
    } else {
      NA_real_
    }
    
    tibble::tibble(
      variable = variable_name,
      r2_linear = r2_linear,
      r2_loess = r2_loess,
      pattern_score = pattern_score
    )
  }
  
  pattern_scores <- df_long |>
    dplyr::group_by(variable) |>
    dplyr::group_split() |>
    purrr::map_dfr(calc_pattern_score)
  
  # Preserve the requested variable order.
  pattern_scores$variable <- factor(
    pattern_scores$variable,
    levels = vars
  )
  
  # ----- Plot trajectories -----
  
  p <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x = .data[[rank_var]],
      y = .data$value
    )
  ) +
    ggplot2::geom_point(
      alpha = point_alpha,
      size = point_size
    ) +
    ggplot2::geom_smooth(
      method = "loess",
      se = show_se,
      color = smooth_color
    ) +
    ggplot2::facet_wrap(
      ~ variable,
      scales = "free_y",
      ncol = ncol
    ) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      x = rank_var,
      y = "Window Score"
    )
  
  # Attach pattern scores to the plot.
  attr(p, "pattern_scores") <- pattern_scores
  
  return(p)
}
