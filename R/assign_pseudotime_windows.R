#' Assign pseudotime windows
#'
#' Assign cells to pseudotime windows based on their pseudotime values,
#' either across all cells or separately within each cell annotation.
#' Window assignments are stored in the Seurat object's metadata.
#'
#' @param se A Seurat object containing the pseudotime values and,
#'   when using stratified mode, the cell annotations.
#' @param pseudotime_col Character string specifying the metadata column
#'   containing numeric pseudotime values.
#' @param mode Character string specifying the window assignment strategy:
#'   "sliding" or "stratified". Defaults to "sliding".
#' @param window_size Positive integer specifying the number of cells per
#'   window in sliding mode. The final window is merged with the preceding
#'   window if it contains fewer cells than the specified window size.
#'   Not used in stratified mode.
#' @param windows_per_anno A named list of positive integers specifying
#'   the number of windows to create for each cell annotation in stratified
#'   mode. Not used in sliding mode.
#' @param cell_anno_col Character string specifying the metadata column
#'   containing cell annotations in stratified mode. Not used in sliding
#'   mode.
#' @param output_col Character string specifying the metadata column in
#'   which to store the window assignments. Defaults to
#'   "pseudotime_windows".
#'
#' @return The Seurat object with cell assignments stored in the metadata
#'   column specified by `output_col`. Cells with missing pseudotime values,
#'   and cells without a matching annotation in stratified mode, receive
#'   `NA`.
#'
#' @details
#' In sliding mode, cells are ordered by pseudotime and assigned to
#' consecutive windows containing up to `window_size` cells. If the final
#' window contains fewer cells than the specified window size, it is
#' merged with the preceding window when more than one window exists.
#'
#' In stratified mode, cells are ordered by pseudotime within each
#' annotation and divided into approximately equal-sized windows.
#' Window identifiers are numbered globally across annotations, following
#' the order of annotations in `windows_per_anno`. The identifiers
#' reserved for an annotation are retained even when that annotation
#' contains no cells with valid pseudotime.
#'
#' Cells are grouped into windows based on their order along pseudotime,
#' with approximately equal numbers of cells per window. The pseudotime
#' range covered by each window may therefore differ.
#'
#' @export
assign_pseudotime_windows <- function(
    se,
    pseudotime_col,
    mode = c("sliding", "stratified"),
    window_size = NULL,
    windows_per_anno = NULL,
    cell_anno_col = NULL,
    output_col = "pseudotime_windows"
) {
  
  # ----- Validate arguments -----
  
  mode <- match.arg(mode)
  
  if (!is.character(pseudotime_col) ||
      length(pseudotime_col) != 1L ||
      is.na(pseudotime_col) ||
      !nzchar(pseudotime_col)) {
    stop("`pseudotime_col` must be a non-empty character string.")
  }
  
  if (!is.character(output_col) ||
      length(output_col) != 1L ||
      is.na(output_col) ||
      !nzchar(output_col)) {
    stop("`output_col` must be a non-empty character string.")
  }
  
  meta <- se@meta.data
  
  if (!pseudotime_col %in% colnames(meta)) {
    stop(
      "The pseudotime column `",
      pseudotime_col,
      "` was not found in the Seurat metadata."
    )
  }
  
  if (!is.numeric(meta[[pseudotime_col]])) {
    stop("`pseudotime_col` must contain numeric values.")
  }
  
  # ----- Initialize output column -----
  
  meta[[output_col]] <- NA_character_
  
  # ----- Sliding mode -----
  
  if (mode == "sliding") {
    
    if (!is.numeric(window_size) ||
        length(window_size) != 1L ||
        is.na(window_size) ||
        !is.finite(window_size) ||
        window_size < 1 ||
        window_size != floor(window_size)) {
      stop(
        "`window_size` must be a positive integer ",
        "when mode = 'sliding'."
      )
    }
    
    # Extract cell identifiers and pseudotime values.
    df <- data.frame(
      cell = rownames(meta),
      pseudotime = meta[[pseudotime_col]],
      stringsAsFactors = FALSE
    )
    
    # Order cells by pseudotime, placing missing values last.
    df_ordered <- df[
      order(df$pseudotime, na.last = TRUE),
      ,
      drop = FALSE
    ]
    
    # Identify cells with valid pseudotime.
    non_na_idx <- which(!is.na(df_ordered$pseudotime))
    n_non_na <- length(non_na_idx)
    
    # Assign windows only when valid pseudotime values are available.
    if (n_non_na > 0L) {
      
      window_ids <- ceiling(
        seq_len(n_non_na) / window_size
      )
      
      # Merge the final window if it is smaller than window_size.
      last_window_id <- max(window_ids)
      
      last_window_size <- sum(
        window_ids == last_window_id
      )
      
      if (last_window_size < window_size &&
          last_window_id > 1L) {
        window_ids[window_ids == last_window_id] <-
          last_window_id - 1L
      }
      
      # Assign window labels.
      labels <- rep(NA_character_, nrow(df_ordered))
      
      labels[non_na_idx] <- paste0(
        "pt_window_",
        window_ids
      )
      
      # Map labels back to the original metadata row order.
      labels_by_cell <- stats::setNames(
        labels,
        df_ordered$cell
      )
      
      meta[[output_col]] <- unname(
        labels_by_cell[rownames(meta)]
      )
    }
    
    # ----- Stratified mode -----
    
  } else if (mode == "stratified") {
    
    if (!is.character(cell_anno_col) ||
        length(cell_anno_col) != 1L ||
        is.na(cell_anno_col) ||
        !nzchar(cell_anno_col)) {
      stop(
        "`cell_anno_col` must be a non-empty character string ",
        "when mode = 'stratified'."
      )
    }
    
    if (!cell_anno_col %in% colnames(meta)) {
      stop(
        "The annotation column `",
        cell_anno_col,
        "` was not found in the Seurat metadata."
      )
    }
    
    if (!is.list(windows_per_anno) ||
        length(windows_per_anno) == 0L ||
        is.null(names(windows_per_anno)) ||
        anyNA(names(windows_per_anno)) ||
        any(!nzchar(names(windows_per_anno))) ||
        anyDuplicated(names(windows_per_anno))) {
      stop(
        "`windows_per_anno` must be a non-empty named list ",
        "with unique, non-empty annotation names."
      )
    }
    
    valid_window_counts <- vapply(
      windows_per_anno,
      function(x) {
        is.numeric(x) &&
          length(x) == 1L &&
          !is.na(x) &&
          is.finite(x) &&
          x >= 1 &&
          x == floor(x)
      },
      logical(1)
    )
    
    if (!all(valid_window_counts)) {
      stop(
        "Each element of `windows_per_anno` must be ",
        "a positive integer."
      )
    }
    
    # Global window counter.
    current_window_id <- 1L
    
    # Process annotations in the supplied order.
    for (anno in names(windows_per_anno)) {
      
      n_windows <- windows_per_anno[[anno]]
      
      # Identify cells in this annotation with valid pseudotime.
      idx <- which(
        !is.na(meta[[cell_anno_col]]) &
          as.character(meta[[cell_anno_col]]) == anno &
          !is.na(meta[[pseudotime_col]])
      )
      
      # Reserve window IDs even if no cells are available.
      if (length(idx) == 0L) {
        current_window_id <- current_window_id + n_windows
        next
      }
      
      # Order cells by pseudotime within the annotation.
      idx_ordered <- idx[
        order(meta[[pseudotime_col]][idx])
      ]
      
      n_cells <- length(idx_ordered)
      
      # Use no more windows than available cells.
      n_assigned_windows <- min(n_windows, n_cells)
      
      # Divide ordered cells into approximately equal-sized groups.
      local_ids <- floor(
        (seq_len(n_cells) - 1) *
          n_assigned_windows / n_cells
      ) + 1L
      
      # Convert local IDs to globally unique window IDs.
      global_ids <- local_ids + current_window_id - 1L
      
      meta[[output_col]][idx_ordered] <- paste0(
        "pt_window_",
        global_ids
      )
      
      # Reserve all requested IDs for this annotation.
      current_window_id <- current_window_id + n_windows
    }
  }
  
  # ----- Return updated Seurat object -----
  
  se@meta.data <- meta
  
  return(se)
}
