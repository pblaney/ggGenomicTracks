#' @name coord_locus
#' @title Cartesian coordinates zoomed to a genomic locus
#'
#' @description
#' Sets the x-axis of a plot to a genomic window. Every layer whose data has a
#' `seqnames` column is subset to rows on the window's chromosome that overlap the
#' window, before aesthetics are evaluated, so discrete scales only show levels present
#' in the window. Layers without a `seqnames` column are left untouched.
#'
#' Also adds an x-axis scale named after the chromosome (e.g. "Chr 3 (Mbp)") with
#' labels in Mbp for windows wider than 10 kbp and bp otherwise. Add your own
#' `scale_x_continuous()` after `coord_locus()` to override it.
#'
#' @param window_of_interest Genomic window as a string, e.g. `"chr3:187734001-187751500"`.
#' @param flank Number of bp to extend both ends of the window, default: 0
#' @param expand If `TRUE`, the default, adds a small expansion to the window limits.
#' @param clip Should drawing be clipped to the extent of the plot panel? Default: `"on"`
#'
#' @returns List of a coord and an x-axis scale, added to a plot with `+`
#' @export
#' @keywords coord
coord_locus <- function(window_of_interest, flank = 0, expand = TRUE, clip = "on") {
  # Parse the window plus flank on both ends
  locus_gr <- gUtils::parse.gr(window_of_interest) + flank
  locus <- list("chrom" = as.character(GenomeInfoDb::seqnames(locus_gr)),
                "start" = GenomicRanges::start(locus_gr),
                "end" = GenomicRanges::end(locus_gr))

  # Build a cartesian coord fixed to the window and extend it with locus subsetting
  base_coord <- ggplot2::coord_cartesian(xlim = c(locus$start, locus$end), expand = expand, clip = clip)
  locus_coord <- ggplot2::ggproto("CoordLocus", base_coord,
    locus = locus,
    setup_data = function(self, data, params = list()) {
      lapply(data, subset_to_locus, locus = self$locus)
    }
  )

  # Pair the coord with an x-axis scale matching the window width
  if((locus$end - locus$start + 1) > 10000) {
    axis_unit <- " (Mbp)"
    axis_labels <- scales::label_comma(accuracy = 0.001, scale = 1e-6)
  } else {
    axis_unit <- " (bp)"
    axis_labels <- scales::label_comma(accuracy = 1)
  }
  axis_name <- paste0(stringr::str_replace(string = locus$chrom, pattern = "^chr", replacement = "Chr "), axis_unit)

  list(locus_coord,
       ggplot2::scale_x_continuous(name = axis_name, labels = axis_labels, expand = c(0.01, 0.01)))
}

#' @name subset_to_locus
#' @title Subset a layer's data to a genomic locus
#'
#' @description
#' Keeps rows on the locus chromosome that overlap the locus. Chromosome names are
#' compared without the `chr` prefix so `1` and `chr1` style data both match. Data
#' without a `seqnames` column is returned unchanged. For internal use only.
#'
#' @param data Data frame of a single layer
#' @param locus List with `chrom`, `start`, and `end` of the locus
#' @keywords internal
subset_to_locus <- function(data, locus) {
  if(!is.data.frame(data) || !"seqnames" %in% names(data)) {
    return(data)
  }

  # Match chromosome regardless of chr prefix
  on_chrom <- stringr::str_remove(as.character(data$seqnames), "^chr") == stringr::str_remove(locus$chrom, "^chr")

  # Match overlap with window, using start alone for single position data
  in_window <- TRUE
  if("start" %in% names(data)) {
    range_end <- if("end" %in% names(data)) data$end else data$start
    in_window <- range_end >= locus$start & data$start <= locus$end
  }

  data[on_chrom & in_window, , drop = FALSE]
}
