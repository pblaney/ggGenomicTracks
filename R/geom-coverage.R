#' @name geom_coverage
#' @title Binned coverage track
#'
#' @description
#' Draws read coverage per genomic bin as points, as in gTrack, e.g. fragCounter or
#' dryclean output. `x` defaults to the middle of each bin and `y` to the
#' `reads.corrected` column. Whole-genome coverage has millions of bins, more than can
#' be seen or drawn quickly, so bins are summarized into `n_bins` bins across the
#' x-axis when there are more than that; zoomed in, every bin is drawn.
#'
#' Also adds a "Coverage" y-axis title, and with `y_range`, a y-axis scale. Add your
#' own `labs()` or scale after it to override them.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()], e.g. `y` for
#'   another coverage column or `colour`.
#' @param data GRanges object or data frame of bins with `seqnames`, `start`, `end`,
#'   and coverage columns.
#' @param n_bins Largest number of points across the x-axis, default: 3000. Use `Inf`
#'   to draw every bin.
#' @param summary Summary of the coverage of bins combined into one point, `"mean"`
#'   (default) or `"median"`.
#' @param y_range Limits of the y-axis, with coverage outside drawn at the limits.
#'   Default: the range of the data.
#' @param point_size Size of the points in mm, default: 0.4
#' @param ... Other arguments passed to [ggplot2::layer()], e.g. `colour` or `alpha`.
#' @param na.rm,show.legend,inherit.aes See [ggplot2::layer()].
#'
#' @returns List of a ggplot2 layer, a y-axis title, and with `y_range`, a y-axis scale,
#'   added to a plot with `+`
#' @export
#' @keywords geom
#'
#' @examples
#' \dontrun{
#' coverage <- readRDS("sample.fragcounter.cov.rds")
#'
#' ggplot2::ggplot(coverage) +
#'   geom_coverage(y_range = c(0, 3)) +
#'   coord_genome() +
#'   theme_genomic()
#' }
geom_coverage <- function(mapping = NULL, data = NULL, n_bins = 3000, summary = c("mean", "median"),
                          y_range = NULL, point_size = 0.4, ...,
                          na.rm = TRUE, show.legend = NA, inherit.aes = TRUE) {
  summary <- match.arg(summary)

  # Default bin middle and coverage aesthetics, overridden by user mapping
  coverage_mapping <- ggplot2::aes(x = (start + end) / 2, y = reads.corrected)
  for(aes_name in names(mapping)) {
    coverage_mapping[[aes_name]] <- mapping[[aes_name]]
  }

  coverage_params <- list(n_bins = n_bins, summary = summary, size = point_size, na.rm = na.rm, ...)
  if(!any(c("colour", "color") %in% c(names(coverage_mapping), names(coverage_params)))) {
    coverage_params$colour <- "gray25"
  }

  coverage_layer <- ggplot2::layer(geom = ggplot2::GeomPoint, stat = StatCoverage,
                                   position = ggplot2::PositionIdentity,
                                   data = as_genomic_df(data), mapping = coverage_mapping,
                                   show.legend = show.legend, inherit.aes = inherit.aes,
                                   params = coverage_params)

  coverage_scale <- NULL
  if(!is.null(y_range)) {
    coverage_scale <- ggplot2::scale_y_continuous(limits = y_range, oob = scales::oob_squish)
  }

  list(coverage_layer, coverage_scale, ggplot2::labs(y = "Coverage"))
}

#' @rdname geom_coverage
#' @format NULL
#' @usage NULL
#' @export
StatCoverage <- ggplot2::ggproto("StatCoverage", ggplot2::Stat,
  required_aes = c("x", "y"),

  compute_panel = function(data, scales, n_bins = 3000, summary = "mean") {
    data <- data[!is.na(data$x) & !is.na(data$y), , drop = FALSE]
    if(nrow(data) <= n_bins) {
      return(data)
    }

    # Combine bins into equal slices of the x-axis, within each group
    x_range <- range(data$x)
    slice <- pmin(floor((data$x - x_range[1]) / diff(x_range) * n_bins), n_bins - 1)
    slice_key <- data$group * n_bins + slice
    first_rows <- !duplicated(slice_key)
    combined <- data[first_rows, , drop = FALSE]
    slice_index <- match(slice_key, slice_key[first_rows])

    counts <- tabulate(slice_index, nbins = nrow(combined))
    combined$x <- as.vector(rowsum(data$x, slice_index, reorder = TRUE)) / counts
    combined$y <- if(summary == "mean") {
      as.vector(rowsum(data$y, slice_index, reorder = TRUE)) / counts
    } else {
      as.vector(tapply(data$y, slice_index, stats::median))
    }
    combined
  }
)
