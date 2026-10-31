#' @name geom_ridge_track
#' @title Signal track
#'
#' @description
#' Draws a quantitative signal over genomic ranges, such as ChIP-seq coverage, a
#' bedGraph, or GISTIC scores, as a shaded area under a step curve. Each range is drawn
#' flat at its value from `start` to `end`. Ranges further apart than `max_gap` are not
#' joined, so missing data is left blank.
#'
#' There are two modes, set by the mapping:
#' * **Overlay**: map the signal to `y` to draw curves on a shared continuous axis,
#'   e.g. `aes(y = score, fill = group)`. The area is shaded between 0 and the signal.
#' * **Ridges**: map a discrete variable to `y` and the signal to `height` to draw one
#'   curve per row, e.g. `aes(y = sample, height = signal)`. Heights are scaled so the
#'   largest absolute value spans `ridge_scale` rows.
#'
#' `xmin` and `xmax` default to the `start` and `end` columns. For single position
#' data, map `x` instead, e.g. `aes(x = position, y = signal)`, and the points are
#' joined directly.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()].
#' @param data GRanges object or data frame of ranges with `seqnames`, `start`, and
#'   `end` columns plus one or more signal columns, in long format with one row per
#'   range and group.
#' @param ridge_scale In ridge mode, the height of the largest absolute signal in
#'   rows, default: 0.9. Values above 1 let ridges overlap.
#' @param max_gap Largest distance in bp between the end of one range and the start of
#'   the next for the two to be joined, default: 1. Use `Inf` to always join.
#' @param outline.type Which sides of the area to outline when `colour` is set, one
#'   of `"upper"`, `"lower"`, `"both"`, or `"full"`. See [ggplot2::geom_ribbon()].
#' @param ... Other arguments passed to [ggplot2::layer()], e.g. `colour`, `fill`,
#'   `alpha`, or `linewidth`.
#' @param na.rm,show.legend,inherit.aes See [ggplot2::layer()].
#'
#' @returns ggplot2 layer
#' @export
#' @keywords geom
#'
#' @examples
#' \dontrun{
#' # GISTIC scores, deletions drawn below zero; segments are ~10 kb apart so join them
#' gistic <- data.table::fread("scores.gistic", check.names = TRUE)
#' gistic <- gistic[, .(seqnames = Chromosome, start = Start, end = End, Type,
#'                      G.score = ifelse(Type == "Del", -G.score, G.score))]
#'
#' ggplot2::ggplot(gistic) +
#'   geom_ridge_track(ggplot2::aes(y = G.score, fill = Type, group = Type), max_gap = 1e5) +
#'   coord_locus("chr13:48303751-48481890", flank = 1e6) +
#'   theme_genomic()
#'
#' # ChIP-seq signal, one ridge per sample
#' chip <- devgru::read_bed_file("signal.bedGraph", additional_col_names = c("sample", "signal"))
#'
#' ggplot2::ggplot(chip) +
#'   geom_ridge_track(ggplot2::aes(y = sample, height = signal, fill = sample)) +
#'   coord_locus("chr1:203200001-203893000") +
#'   theme_genomic()
#' }
geom_ridge_track <- function(mapping = NULL, data = NULL, ridge_scale = 0.9, max_gap = 1,
                             outline.type = "upper", ...,
                             na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  # Default range aesthetics unless single positions are mapped to x, overridden by
  # user mapping
  ridge_mapping <- if("x" %in% names(mapping)) ggplot2::aes() else ggplot2::aes(xmin = start, xmax = end)
  for(aes_name in names(mapping)) {
    ridge_mapping[[aes_name]] <- mapping[[aes_name]]
  }

  ggplot2::layer(geom = GeomRidgeTrack, stat = StatRidgeTrack, position = ggplot2::PositionIdentity,
                 data = as_genomic_df(data), mapping = ridge_mapping,
                 show.legend = show.legend, inherit.aes = inherit.aes,
                 params = list(ridge_scale = ridge_scale, max_gap = max_gap,
                               outline.type = outline.type, na.rm = na.rm, ...))
}

#' @rdname geom_ridge_track
#' @format NULL
#' @usage NULL
#' @export
StatRidgeTrack <- ggplot2::ggproto("StatRidgeTrack", ggplot2::Stat,
  required_aes = "y",
  optional_aes = c("x", "xmin", "xmax", "height"),
  dropped_aes = c("xmin", "xmax", "height"),

  setup_params = function(data, params) {
    # Ridge heights are scaled by the largest absolute signal across the whole layer
    # so all rows share one scale
    params$height_max <- if("height" %in% names(data)) max(abs(data$height), na.rm = TRUE) else NA
    params
  },

  compute_group = function(data, scales, ridge_scale = 0.9, max_gap = 1, height_max = NA) {
    # Single positions are joined directly
    if(!"xmin" %in% names(data)) {
      data$xmin <- data$x
      data$xmax <- data$x
      max_gap <- Inf
    }
    data <- data[order(data$xmin), , drop = FALSE]

    # Baseline and top of the area, in rows for ridges or signal units for overlay
    if("height" %in% names(data)) {
      baseline <- data$y
      top <- data$y + data$height / height_max * ridge_scale
    } else {
      baseline <- rep(0, nrow(data))
      top <- data$y
    }

    # Each range becomes a flat step, held up to the start of the next range when the
    # two are joined, with an NA row after each range followed by a gap so the area
    # breaks there
    n <- nrow(data)
    next_start <- c(data$xmin[-1], NA)
    gap_after <- c(data$xmin[-1] - data$xmax[-n] > max_gap, TRUE)
    step_end <- ifelse(gap_after | is.na(next_start), data$xmax, next_start)
    gap_after[n] <- FALSE
    row_index <- rep(seq_len(n), each = 2)
    steps <- data[row_index, setdiff(names(data), c("x", "xmin", "xmax", "y", "height")), drop = FALSE]
    steps$x <- as.vector(rbind(data$xmin, step_end))
    steps$ymin <- baseline[row_index]
    steps$ymax <- top[row_index]
    steps$y <- steps$ymax
    steps$break_after <- rep(c(FALSE, TRUE), n) & gap_after[row_index]

    break_rows <- steps[steps$break_after, , drop = FALSE]
    if(nrow(break_rows) > 0) {
      break_rows$ymin <- NA
      break_rows$ymax <- NA
      break_rows$y <- NA
      steps <- rbind(steps, break_rows)
      steps <- steps[order(steps$x, is.na(steps$ymax)), , drop = FALSE]
    }
    steps$break_after <- NULL
    steps
  }
)

#' @rdname geom_ridge_track
#' @format NULL
#' @usage NULL
#' @export
GeomRidgeTrack <- ggplot2::ggproto("GeomRidgeTrack", ggplot2::GeomRibbon,
  default_aes = ggplot2::aes(colour = NA, fill = "darkred", linewidth = 0.3,
                             linetype = 1, alpha = NA),

  # NA rows from the stat mark gaps between ranges; split the area there without the
  # missing value warning
  handle_na = function(self, data, params) {
    missing <- !stats::complete.cases(data[c("x", "ymin", "ymax")])
    if(!any(missing)) {
      return(data)
    }
    data$group <- cumsum(c(TRUE, missing[-1] != missing[-length(missing)] |
                                 data$group[-1] != data$group[-nrow(data)]))
    data[!missing, , drop = FALSE]
  }
)
