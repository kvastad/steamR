#' Order Fisher test results
#'
#' Reorders Fisher test results according to a user-defined sequence
#' and converts the grouping column to an ordered factor.
#'
#' @param fisher_results Data frame containing Fisher test results.
#' @param order_levels Character vector specifying the desired order
#'   of the groups in the `cluster` column.
#' @param group_col Name of the column containing the group labels.
#'   Defaults to "cluster".
#'
#' @return A data frame with rows ordered according to `order_levels`.
#'
#' @export

order_fisher_results <- function(
    fisher_results,
    order_levels,
    group_col = "cluster"
) {
  
  # -------------------------------
  # Check inputs
  # -------------------------------
  
  if (!is.data.frame(fisher_results)) {
    stop("`fisher_results` must be a data frame.")
  }
  
  if (!group_col %in% colnames(fisher_results)) {
    stop(
      "The grouping column `",
      group_col,
      "` was not found in `fisher_results`."
    )
  }
  
  if (!is.character(order_levels) || length(order_levels) == 0) {
    stop("`order_levels` must be a non-empty character vector.")
  }
  
  if (anyNA(order_levels) || anyDuplicated(order_levels)) {
    stop("`order_levels` must not contain missing or duplicated values.")
  }
  
  # -------------------------------
  # Check group labels
  # -------------------------------
  
  group_values <- as.character(fisher_results[[group_col]])
  
  if (anyNA(group_values)) {
    stop("The grouping column contains missing values.")
  }
  
  unexpected_groups <- setdiff(
    unique(group_values),
    order_levels
  )
  
  if (length(unexpected_groups) > 0) {
    stop(
      "The following groups are missing from `order_levels`: ",
      paste(unexpected_groups, collapse = ", ")
    )
  }
  
  # -------------------------------
  # Reorder results
  # -------------------------------
  
  fisher_results[[group_col]] <- factor(
    group_values,
    levels = order_levels,
    ordered = TRUE
  )
  
  fisher_results <- fisher_results[
    order(fisher_results[[group_col]]),
    ,
    drop = FALSE
  ]
  
  return(fisher_results)
}