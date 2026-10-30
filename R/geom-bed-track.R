#' @name geom_bed_track
#' @title Genomic interval track
#'
#' @description
#' Draws BED-like intervals as bars of fixed height, one row per level of the `y`
#' aesthetic. `xmin` and `xmax` default to the `start` and `end` columns. Pair with
#' [coord_locus()] to zoom to a window and [theme_genomic()] for track styling.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()]. Typically
#'   `y` (e.g. `sample`) and optionally `fill`.
#' @param data GRanges object or data frame of intervals with `seqnames`, `start`, and
#'   `end` columns, e.g. from `devgru::read_bed_file()`.
#' @param bar_height Height of the bars in mm, default: 5
#' @param ... Other arguments passed to [gggenes::geom_gene_arrow()], e.g. `colour`
#'   or `alpha`. Bars have no outline unless `colour` is mapped or set.
#' @param na.rm,show.legend,inherit.aes See [ggplot2::layer()].
#'
#' @returns ggplot2 layer
#' @export
#' @keywords geom
#'
#' @examples
#' \dontrun{
#' chrom_states <- devgru::read_bed_file("segments.bed", additional_col_names = "segment_label")
#'
#' ggplot2::ggplot() +
#'   geom_bed_track(data = chrom_states, ggplot2::aes(y = sample, fill = segment_label)) +
#'   coord_locus("chr1:203296001-203313000") +
#'   theme_genomic()
#'
#' # Tornado sort: order rows by interval width, largest at top
#' ggplot2::ggplot() +
#'   geom_bed_track(data = chrom_states, ggplot2::aes(y = reorder(sample, end - start))) +
#'   coord_locus("chr1:203296001-203313000")
#' }
geom_bed_track <- function(mapping = NULL, data = NULL, bar_height = 5, ...,
                           na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  # Default interval aesthetics, overridden by user mapping
  bed_mapping <- ggplot2::aes(xmin = start, xmax = end)
  for(aes_name in names(mapping)) {
    bed_mapping[[aes_name]] <- mapping[[aes_name]]
  }

  # No outline by default so adjacent intervals read as one bar
  dots <- list(...)
  if(!any(c("colour", "color") %in% c(names(bed_mapping), names(dots)))) {
    dots$colour <- NA
  }

  do.call(gggenes::geom_gene_arrow,
          c(list(mapping = bed_mapping,
                 data = as_genomic_df(data),
                 arrow_body_height = ggplot2::unit(bar_height, "mm"),
                 arrowhead_height = ggplot2::unit(bar_height, "mm"),
                 arrowhead_width = ggplot2::unit(0, "mm"),
                 na.rm = na.rm,
                 show.legend = show.legend,
                 inherit.aes = inherit.aes),
            dots))
}

#' @name as_genomic_df
#' @title Convert genomic data to a plain data frame for a ggplot2 layer
#'
#' @description
#' Converts GRanges objects to a data frame with `seqnames`, `start`, `end`, `strand`,
#' and `width` columns plus metadata columns. Data frames and `NULL` are returned
#' unchanged. For internal use only.
#'
#' @param data GRanges object, data frame, or `NULL`
#' @keywords internal
as_genomic_df <- function(data) {
  if(inherits(data, "GRanges")) {
    data <- as.data.frame(gUtils::gr2dt(data))
  }
  data
}
