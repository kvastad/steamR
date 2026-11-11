#' Plot cell-type composition along a ranked trajectory
#'
#' Visualizes cell-type composition across cell ranks, with pseudotime
#' window boundaries and labels. Optionally displays Fisher's exact
#' test results as odds-ratio annotations and significance stars.
#'
#' Cell types are detected automatically from `cell_anno_col` when
#' `composition_cols` is `NULL`. If `composition_cols` is supplied,
#' cell-type composition is calculated from the corresponding count
#' columns in the Seurat metadata (for example, MetaCell data).
#'
#' @param se Seurat object containing cell annotations or cell-type
#'   composition counts, cell ranks, and pseudotime window assignments.
#' @param fisher_results Optional data frame containing `cluster`,
#'   `adj_p_value`, and `odds_ratio`.
#' @param cell_colors Optional named vector of colors for each cell type.
#'   If `NULL`, ggplot2's default discrete colors are used.
#' @param cell_anno_col Metadata column containing cell-type annotations.
#'   Used when `composition_cols` is `NULL`.
#' @param composition_cols Optional character vector specifying metadata
#'   columns containing cell-type counts. Use this for MetaCell data.
#' @param rank_col Metadata column containing ordered cell ranks.
#' @param window_col Metadata column containing pseudotime window assignments.
#' @param window_label Display window labels as `"number"` or `"full"`.
#' @param title Plot title.
#'
#' @return A ggplot object.
#'
#' @importFrom dplyr %>%
#' @export

plot_cell_composition <- function(
    se,
    fisher_results = NULL,
    cell_colors = NULL,
    cell_anno_col = "cell_anno",
    composition_cols = NULL,
    rank_col = "cell_rank_ordered",
    window_col = "pseudotime_windows",
    window_label = c("number", "full"),
    title = NULL
) {

  # -------------------------------
  # Check inputs
  # -------------------------------

  window_label <- match.arg(window_label)

  meta <- se[[]]

  required_cols <- c(rank_col, window_col)

  if (is.null(composition_cols)) {
    required_cols <- c(required_cols, cell_anno_col)
  } else {
    if (!is.character(composition_cols) || length(composition_cols) == 0) {
      stop("`composition_cols` must be a non-empty character vector.")
    }

    if (anyNA(composition_cols) || anyDuplicated(composition_cols)) {
      stop("`composition_cols` must not contain missing or duplicated values.")
    }

    required_cols <- c(required_cols, composition_cols)
  }

  if (!all(required_cols %in% colnames(meta))) {
    stop(
      "Missing required metadata columns: ",
      paste(setdiff(required_cols, colnames(meta)), collapse = ", ")
    )
  }

  # -------------------------------
  # Detect cell types
  # -------------------------------

  if (is.null(composition_cols)) {

    cell_anno <- meta[[cell_anno_col]]

    # Preserve factor level order when available.
    # For other data types, use the order of first appearance.
    if (is.factor(cell_anno)) {

      cell_types <- levels(cell_anno)
      cell_types <- cell_types[
        cell_types %in% as.character(cell_anno)
      ]

    } else {

      cell_types <- unique(as.character(cell_anno))
      cell_types <- cell_types[!is.na(cell_types)]
    }

  } else {

    cell_types <- composition_cols
  }

  if (length(cell_types) == 0) {
    stop("No cell types found.")
  }

  # Check user-provided cell colors, if supplied.
  if (!is.null(cell_colors)) {

    if (
      is.null(names(cell_colors)) ||
      any(is.na(names(cell_colors))) ||
      any(names(cell_colors) == "") ||
      anyDuplicated(names(cell_colors))
    ) {
      stop("`cell_colors` must be a uniquely named color vector.")
    }

    if (!all(cell_types %in% names(cell_colors))) {
      stop(
        "Each cell type must have a corresponding color ",
        "in `cell_colors`."
      )
    }

    # Keep only colors corresponding to cell types in the data,
    # and arrange them in the same order as the cell types.
    cell_colors <- cell_colors[cell_types]
  }

  # -------------------------------
  # Check MetaCell composition data
  # -------------------------------

  if (!is.null(composition_cols)) {

    composition_data <- meta[, composition_cols, drop = FALSE]

    if (!all(vapply(composition_data, is.numeric, logical(1)))) {
      stop("All columns specified in `composition_cols` must be numeric.")
    }

    composition_matrix <- as.matrix(composition_data)

    if (anyNA(composition_matrix) || any(!is.finite(composition_matrix))) {
      stop("Columns specified in `composition_cols` must not contain missing or non-finite values.")
    }

    if (any(composition_matrix < 0)) {
      stop("Cell-type counts in `composition_cols` must not be negative.")
    }
  }

  # -------------------------------
  # Check Fisher results
  # -------------------------------

  if (!is.null(fisher_results)) {

    required_fisher_cols <- c("cluster", "adj_p_value", "odds_ratio")

    if (!all(required_fisher_cols %in% colnames(fisher_results))) {
      stop(
        "`fisher_results` must contain: ",
        paste(required_fisher_cols, collapse = ", ")
      )
    }

    if (any(
      !is.na(fisher_results$odds_ratio) &
        fisher_results$odds_ratio <= 0
    )) {
      stop("All non-missing odds ratios must be greater than zero.")
    }
  }

  # -------------------------------
  # Prepare cell ranks
  # -------------------------------

  rank_values <- meta[[rank_col]]

  # Preserve existing factor ordering.
  if (is.factor(rank_values)) {

    rank_levels <- levels(rank_values)
    rank_levels <- rank_levels[
      rank_levels %in% as.character(rank_values)
    ]

  } else if (is.numeric(rank_values)) {

    rank_levels <- as.character(sort(unique(
      rank_values[!is.na(rank_values)]
    )))

  } else {

    rank_levels <- unique(as.character(rank_values))
    rank_levels <- rank_levels[!is.na(rank_levels)]

    # Sort numerically if ranks are numeric strings.
    numeric_ranks <- suppressWarnings(as.numeric(rank_levels))

    if (length(numeric_ranks) > 0 && all(!is.na(numeric_ranks))) {
      rank_levels <- rank_levels[order(numeric_ranks)]
    }
  }

  if (length(rank_levels) == 0) {
    stop("No valid cell ranks found.")
  }

  # -------------------------------
  # Prepare pseudotime windows
  # -------------------------------

  window_values <- unique(
    as.character(meta[[window_col]][!is.na(meta[[window_col]])])
  )

  if (length(window_values) == 0) {
    stop("No pseudotime windows found.")
  }

  # Order windows numerically when their names end in a number
  # (for example, window_1, window_2, window_10).
  window_numbers <- suppressWarnings(
    as.numeric(sub(".*_", "", window_values))
  )

  if (all(!is.na(window_numbers))) {

    window_values <- window_values[order(window_numbers)]

  } else if (is.factor(meta[[window_col]])) {

    window_levels <- levels(meta[[window_col]])

    window_values <- window_levels[
      window_levels %in% window_values
    ]
  }

  # Create internal plotting metadata.
  meta_plot <- data.frame(
    .rank_label = as.character(rank_values),
    .rank = match(as.character(rank_values), rank_levels),
    .window = as.character(meta[[window_col]])
  )

  meta_plot$.window <- factor(
    meta_plot$.window,
    levels = window_values,
    ordered = TRUE
  )

  # -------------------------------
  # Cell-type composition
  # -------------------------------

  if (is.null(composition_cols)) {

    # Single-cell data: count cells in each cell type per rank.
    meta_plot$.celltype <- as.character(cell_anno)

    plot_long <- meta_plot %>%
      dplyr::filter(
        .celltype %in% cell_types,
        !is.na(.rank)
      ) %>%
      dplyr::count(
        .rank,
        .celltype,
        name = "Count"
      ) %>%
      dplyr::group_by(.rank) %>%
      dplyr::mutate(
        Fraction = Count / sum(Count)
      ) %>%
      dplyr::ungroup()

  } else {

    # MetaCell data: use the supplied cell-type count columns.
    # Calculate each cell type's fraction within each MetaCell.
    composition_totals <- rowSums(composition_matrix)

    plot_long <- data.frame(
      .rank = rep(meta_plot$.rank, times = length(cell_types)),
      .celltype = rep(cell_types, each = nrow(meta_plot)),
      Count = as.vector(composition_matrix),
      Total = rep(composition_totals, times = length(cell_types))
    ) %>%
      dplyr::filter(
        !is.na(.rank),
        Total > 0
      ) %>%
      dplyr::mutate(
        Fraction = Count / Total
      )
  }

  if (nrow(plot_long) == 0) {
    stop("No valid cell-type composition data and cell ranks were found.")
  }

  plot_long$CellType <- factor(
    plot_long$.celltype,
    levels = cell_types
  )

  # -------------------------------
  # Pseudotime window boundaries
  # -------------------------------

  window_bounds <- meta_plot %>%
    dplyr::filter(
      !is.na(.window),
      !is.na(.rank)
    ) %>%
    dplyr::group_by(.window) %>%
    dplyr::summarise(
      start = min(.rank),
      end = max(.rank),
      .groups = "drop"
    ) %>%
    dplyr::arrange(start)

  window_bounds$pseudotime_windows <- as.character(
    window_bounds$.window
  )

  window_labels <- window_bounds %>%
    dplyr::mutate(
      xpos = (start + end) / 2,
      label = if (window_label == "number") {
        sub(".*_", "", pseudotime_windows)
      } else {
        pseudotime_windows
      }
    )

  # Dashed lines between windows.
  vline_positions <- if (nrow(window_bounds) > 1) {
    window_bounds$end[-nrow(window_bounds)] + 0.5
  } else {
    numeric(0)
  }

  # X-axis labels: last rank of each window.
  x_breaks <- window_bounds$end
  x_labels <- rank_levels[x_breaks]

  # -------------------------------
  # Base plot
  # -------------------------------

  has_fisher <- !is.null(fisher_results)

  window_y <- if (has_fisher) 1.17 else 1.08
  y_max <- if (has_fisher) 1.32 else 1.20

  p <- ggplot2::ggplot(
    plot_long,
    ggplot2::aes(
      x = .rank,
      y = Fraction,
      fill = CellType
    )
  ) +

    # Cell composition bars.
    ggplot2::geom_bar(
      stat = "identity",
      width = 1
    ) +

    # Pseudotime window boundaries.
    ggplot2::geom_vline(
      xintercept = vline_positions,
      linetype = "dashed",
      color = "black"
    ) +

    # Pseudotime window labels.
    ggplot2::geom_text(
      data = window_labels,
      ggplot2::aes(
        x = xpos,
        y = window_y,
        label = label
      ),
      inherit.aes = FALSE,
      angle = if (window_label == "full") 90 else 0,
      hjust = if (window_label == "full") 0 else 0.5,
      vjust = if (window_label == "full") 0 else -0.3,
      size = if (window_label == "full") 4 else 5
    )

  # Use manual cell-type colors only when supplied.
  if (!is.null(cell_colors)) {

    p <- p +
      ggplot2::scale_fill_manual(
        values = cell_colors,
        name = "Cell Type",
        drop = FALSE
      )

  } else {

    p <- p +
      ggplot2::scale_fill_discrete(
        name = "Cell Type",
        drop = FALSE
      )
  }

  p <- p +

    # X-axis.
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = c(0, 0),
      limits = c(0.5, max(plot_long$.rank) + 0.5)
    ) +

    # Y-axis.
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.5),
      expand = c(0, 0)
    ) +

    ggplot2::coord_cartesian(
      ylim = c(0, y_max),
      clip = "off"
    ) +

    ggplot2::labs(
      x = if (is.null(composition_cols)) "Cell rank" else "MetaCell rank",
      y = "Fraction of cells",
      title = title
    ) +

    ggplot2::theme_classic() +

    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = 15
      ),
      axis.text.y = ggplot2::element_text(size = 15),
      axis.title.x = ggplot2::element_text(size = 16),
      axis.title.y = ggplot2::element_text(size = 16),
      legend.title = ggplot2::element_text(size = 15),
      legend.text = ggplot2::element_text(size = 15),
      plot.title = ggplot2::element_text(
        size = 20,
        face = "bold",
        margin = ggplot2::margin(b = 25)
      ),
      plot.margin = ggplot2::margin(
        t = 25, r = 10, b = 10, l = 10
      )
    )

  # -------------------------------
  # Optional: Fisher annotations
  # -------------------------------

  if (has_fisher) {

    fisher_plot <- fisher_results %>%
      dplyr::select(
        cluster,
        adj_p_value,
        odds_ratio
      ) %>%
      dplyr::inner_join(
        window_bounds %>%
          dplyr::select(
            pseudotime_windows,
            start,
            end
          ),
        by = c("cluster" = "pseudotime_windows")
      ) %>%
      dplyr::mutate(
        xmin = start - 0.4,
        xmax = end + 0.4,
        xmid = (xmin + xmax) / 2,
        log2_OR = log2(odds_ratio),
        stars = dplyr::case_when(
          is.na(adj_p_value) ~ "",
          adj_p_value < 0.0001 ~ "****",
          adj_p_value < 0.001 ~ "***",
          adj_p_value < 0.01 ~ "**",
          adj_p_value < 0.05 ~ "*",
          TRUE ~ ""
        )
      )

    # Odds-ratio color annotations.
    p <- p +
      ggnewscale::new_scale_fill() +

      ggplot2::geom_rect(
        data = fisher_plot,
        ggplot2::aes(
          xmin = xmin,
          xmax = xmax,
          ymin = 1.02,
          ymax = 1.07,
          fill = log2_OR
        ),
        inherit.aes = FALSE,
        color = NA
      ) +

      ggplot2::scale_fill_viridis_c(
        option = "D",
        direction = 1,
        name = "log2(odds_ratio)"
      ) +

      # Significance stars.
      ggplot2::geom_text(
        data = fisher_plot,
        ggplot2::aes(
          x = xmid,
          y = 1.10,
          label = stars
        ),
        inherit.aes = FALSE,
        size = 8
      )
  }

  return(p)
}
