#' @keywords internal
"_PACKAGE"

#' @import data.table
#' @import devgru
#' @import ggplot2
#' @import scales
#' @import stringr
#' @importFrom dplyr %>% .data all_of
#' @importFrom gUtils %Q% %&% gr2dt dt2gr
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 Hsapiens
#' @importFrom gggenes geom_gene_arrow
NULL

# Column names referenced via non-standard evaluation in dplyr/ggplot2 calls
utils::globalVariables(c("seqnames", "start", "end", "width", "n",
                         "gene_name", "segment", "orientation", "label", "y_axis", "side_of_anno",
                         "VAF", "bits", "chromhmm_state", "thickness", "sv_type", "colors",
                         "x", "y", "arm", "gieStain",
                         "forward", "arrow_from", "arrow_to", "label_x", "label_hjust", "count", "score", "state"))
