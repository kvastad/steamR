#' Internal helper for Seurat v4/v5 data access
#'
#' Uses `layer` for SeuratObject >= 5, else `slot`. Always prefers "RNA" if it exists,
#' otherwise falls back to "Spatial" (for v5) or the default assay (for v4).
#' Prints which assay is being used.
#' @keywords internal
get_assay_data <- function(se, what = c("counts", "data", "scale.data")) {
  what <- match.arg(what)
  seurat_ver <- utils::packageVersion("SeuratObject")
  
  # Determine assay
  if (seurat_ver >= "5.0.0") {
    assays_available <- SeuratObject::Assays(se)
    if ("Spatial" %in% assays_available) {
      assay <- "Spatial"
    } else if ("RNA" %in% assays_available) {
      assay <- "RNA"
    } else {
      assay <- SeuratObject::DefaultAssay(se)
    }
  } else {
    assays_available <- SeuratObject::Assays(se)
    if ("RNA" %in% assays_available) {
      assay <- "RNA"
    } else {
      assay <- SeuratObject::DefaultAssay(se)
    }
  }
  
  message("Using assay: ", assay, " | slot/layer: ", what)
  
  # Access data
  if (seurat_ver >= "5.0.0") {
    SeuratObject::GetAssayData(se, assay = assay, layer = what)
  } else {
    SeuratObject::GetAssayData(se, assay = assay, slot = what)
  }
}
