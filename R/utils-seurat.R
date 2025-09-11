#' Internal helper for Seurat v4/v5 data access
#'
#' Uses `layer` for SeuratObject >= 5, else `slot`.
#' @keywords internal
get_assay_data <- function(se, assay = SeuratObject::DefaultAssay(se), what = c("counts","data","scale.data")) {
  what <- match.arg(what)
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    SeuratObject::GetAssayData(se, assay = assay, layer = what)
  } else {
    SeuratObject::GetAssayData(se, assay = assay, slot  = what)
  }
}


