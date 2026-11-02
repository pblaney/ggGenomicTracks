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
#' `coord_locus()` is [coord_loci()] with a single window, and x positions are the
#' genomic coordinates themselves.
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
  coord_loci(windows = window_of_interest, flank = flank, expand = expand, clip = clip)
}

#' @name coord_loci
#' @title Cartesian coordinates across multiple genomic loci
#'
#' @description
#' Places several genomic windows side by side on one linear x-axis, so tracks can show
#' distant loci, on the same or different chromosomes, in a single panel. Every layer
#' whose data has `seqnames` and `start` columns has its positions converted to the
#' linear axis before aesthetics are evaluated: rows outside all windows are dropped,
#' and with more than one window, ranges are clipped to the window edges and split
#' across windows they span. Data with two anchors, with `seqnames1`, `start1`,
#' `seqnames2`, and `start2` columns as from `devgru::read_bedpe_file(output_type =
#' "dt")`, has both anchors converted and keeps rows with at least one anchor in view;
#' the anchor out of view is set to `NA`. Layers without these columns are left
#' untouched.
#'
#' Also adds an x-axis scale with breaks and labels in each window's own coordinates,
#' and window labels on the top axis.
#'
#' @param windows Genomic windows as one of:
#'   * A character vector of windows, e.g. `"chr8:127600001-128000000"`. Names, if
#'     given, are used as window labels.
#'   * A data frame with `seqnames`, `start`, and `end` columns and an optional `label`
#'     column.
#'   * A GRanges object, with an optional `label` metadata column.
#' @param flank Number of bp to extend both ends of each window, default: 0
#' @param widths How the axis is shared between windows: `"proportional"` (default) to
#'   each window's length in bp, or `"equal"` so small windows stay readable next to
#'   large ones, each then at its own scale.
#' @param gap Space between windows as a fraction of the total axis, default: 0.02
#' @param dividers Draw dashed lines between windows, default: TRUE
#' @param window_labels Labels of the windows on the top axis, as one of:
#'   * `NULL` (default): labeled with more than one window, by the names of `windows`
#'     or else the chromosome.
#'   * A character vector with one label per window.
#'   * `TRUE` or `FALSE` to show or hide the default labels. In a stack, keep the
#'     labels on the top track only.
#' @inheritParams coord_locus
#'
#' @returns List of a coord and an x-axis scale, added to a plot with `+`
#' @export
#' @keywords coord
#'
#' @examples
#' \dontrun{
#' windows <- c(MYC = "chr8:127700001-128000000", IGH = "chr14:105550001-106880000")
#'
#' ggplot2::ggplot() +
#'   geom_gene_track(genes_of_interest = c("MYC", "PVT1")) +
#'   coord_loci(windows) +
#'   theme_genomic(y_axis = "blank")
#' }
coord_loci <- function(windows, flank = 0, widths = c("proportional", "equal"), gap = 0.02,
                       dividers = TRUE, window_labels = NULL, expand = TRUE, clip = "on") {
  widths <- match.arg(widths)
  layout <- locus_layout(as_locus_windows(windows, flank = flank), widths = widths, gap = gap)

  # Window labels given directly replace the default names
  if(is.character(window_labels)) {
    if(length(window_labels) != nrow(layout)) {
      cli::cli_abort("{.arg window_labels} must have one label per window ({nrow(layout)}), not {length(window_labels)}.")
    }
    layout$label <- window_labels
    window_labels <- TRUE
  }
  if(is.null(window_labels)) {
    window_labels <- nrow(layout) > 1
  }

  # Axis breaks within each window, labeled in that window's coordinates
  breaks <- locus_breaks(layout)
  axis_unit <- if(all(layout$end - layout$start + 1 > 10000)) " (Mbp)" else " (bp)"
  axis_name <- if(nrow(layout) == 1) {
    paste0(stringr::str_replace(layout$chrom, "^chr", "Chr "), axis_unit)
  } else {
    paste0("Position", axis_unit)
  }
  window_axis <- if(window_labels) {
    ggplot2::sec_axis(~ ., breaks = (layout$lin_start + layout$lin_end) / 2, labels = layout$label)
  } else {
    ggplot2::waiver()
  }

  spec <- list(fun = "coord_loci",
               args = list(windows = layout[, c("chrom", "start", "end", "label")], widths = widths, gap = gap,
                           dividers = dividers, expand = expand, clip = clip))
  names(spec$args$windows)[1] <- "seqnames"

  list(new_coord_loci(layout, expand = expand, clip = clip, dividers = dividers, spec = spec),
       ggplot2::scale_x_continuous(name = axis_name, breaks = breaks$lin, labels = breaks$label,
                                   expand = c(0.01, 0.01), sec.axis = window_axis))
}

#' @name coord_genome
#' @title Cartesian coordinates across the whole genome
#'
#' @description
#' Places whole chromosomes end to end on one linear x-axis, for genome-wide tracks
#' such as copy number. Works like [coord_loci()] with every chromosome as a window:
#' the x-axis is labeled with chromosome names, and alternate chromosomes are shaded.
#'
#' @param chromosomes Chromosomes to show, in order, with or without the `chr` prefix.
#'   Default: 1-22, X, and Y.
#' @param seqlengths Named vector of chromosome lengths, default: hg38 from
#'   BSgenome.Hsapiens.UCSC.hg38.
#' @param gap Space between chromosomes as a fraction of the total axis, default: 0
#' @param dividers Draw dashed lines between chromosomes, default: TRUE
#' @param shade Fill color of every other chromosome, default: `"gray95"`. Use `NA` for
#'   none.
#' @inheritParams coord_locus
#'
#' @returns List of a coord and an x-axis scale, added to a plot with `+`
#' @export
#' @keywords coord
#'
#' @examples
#' \dontrun{
#' ggplot2::ggplot(cn_segments) +
#'   geom_bed_track(ggplot2::aes(y = sample, fill = cn)) +
#'   coord_genome() +
#'   theme_genomic()
#' }
coord_genome <- function(chromosomes = c(1:22, "X", "Y"), seqlengths = NULL, gap = 0,
                         dividers = TRUE, shade = "gray95", expand = FALSE, clip = "on") {
  if(is.null(seqlengths)) {
    seqlengths <- GenomeInfoDb::seqlengths(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)
  }
  chrom_names <- paste0("chr", stringr::str_remove(as.character(chromosomes), "^chr"))
  length_names <- paste0("chr", stringr::str_remove(names(seqlengths), "^chr"))
  chrom_lengths <- seqlengths[match(chrom_names, length_names)]
  if(anyNA(chrom_lengths)) {
    cli::cli_abort("No length found for chromosome{?s} {.val {chromosomes[is.na(chrom_lengths)]}} in {.arg seqlengths}.")
  }

  windows <- data.frame(seqnames = chrom_names, start = 1, end = as.numeric(chrom_lengths),
                        label = stringr::str_remove(chrom_names, "^chr"))
  layout <- locus_layout(windows, widths = "proportional", gap = gap)

  spec <- list(fun = "coord_genome",
               args = list(chromosomes = chromosomes, seqlengths = seqlengths, gap = gap, dividers = dividers,
                           shade = shade, expand = expand, clip = clip))

  list(new_coord_loci(layout, expand = expand, clip = clip, dividers = dividers, shade = shade, spec = spec),
       ggplot2::scale_x_continuous(name = "Chromosome", breaks = (layout$lin_start + layout$lin_end) / 2,
                                   labels = layout$label, expand = c(0, 0)))
}

#' @name new_coord_loci
#' @title Build the multi-locus coord
#'
#' @description
#' Extends [ggplot2::coord_cartesian()] with a window layout used to convert layer data
#' to the linear axis, and draws window dividers and shading. For internal use only.
#'
#' @param layout Window layout from [locus_layout()]
#' @param dividers Draw dashed lines between windows
#' @param shade Fill color of every other window, or `NA` for none
#' @param spec List of the coord function name and its arguments, used by
#'   [stack_tracks()] to rebuild the coord for each track
#' @inheritParams coord_locus
#' @keywords internal
new_coord_loci <- function(layout, expand = TRUE, clip = "on", dividers = FALSE, shade = NA, spec = NULL) {
  base_coord <- ggplot2::coord_cartesian(xlim = c(min(layout$lin_start), max(layout$lin_end)),
                                         expand = expand, clip = clip)
  ggplot2::ggproto("CoordLoci", base_coord,
    layout = layout,
    dividers = dividers,
    shade = shade,
    spec = spec,

    setup_data = function(self, data, params = list()) {
      lapply(data, linearize_layer, layout = self$layout)
    },

    # Shading of alternate windows behind the panel grid
    render_bg = function(self, panel_params, theme) {
      background <- ggplot2::ggproto_parent(ggplot2::CoordCartesian, self)$render_bg(panel_params, theme)
      if(is.na(self$shade) || nrow(self$layout) < 2) {
        return(background)
      }
      shaded <- self$layout[seq(2, nrow(self$layout), by = 2), , drop = FALSE]
      x_range <- panel_params$x.range
      grid::grobTree(background,
                     grid::rectGrob(x = grid::unit((shaded$lin_start - x_range[1]) / diff(x_range), "npc"),
                                    width = grid::unit((shaded$lin_end - shaded$lin_start) / diff(x_range), "npc"),
                                    just = "left", gp = grid::gpar(fill = self$shade, col = NA)))
    },

    # Dashed lines between adjacent windows, centered in any gap, over the data
    render_fg = function(self, panel_params, theme) {
      foreground <- ggplot2::ggproto_parent(ggplot2::CoordCartesian, self)$render_fg(panel_params, theme)
      if(!self$dividers || nrow(self$layout) < 2) {
        return(foreground)
      }
      x_range <- panel_params$x.range
      boundaries <- (self$layout$lin_end[-nrow(self$layout)] + self$layout$lin_start[-1]) / 2
      grid::grobTree(foreground,
                     grid::segmentsGrob(x0 = grid::unit((boundaries - x_range[1]) / diff(x_range), "npc"),
                                        x1 = grid::unit((boundaries - x_range[1]) / diff(x_range), "npc"),
                                        y0 = grid::unit(0, "npc"), y1 = grid::unit(1, "npc"),
                                        gp = grid::gpar(col = "gray50", lwd = 0.6, lty = "dashed")))
    }
  )
}

#' @name as_locus_windows
#' @title Standardize genomic windows
#'
#' @description
#' Converts window strings, data frames, or GRanges to a data frame with `seqnames`,
#' `start`, `end`, and `label` columns, extended by the flank. For internal use only.
#'
#' @inheritParams coord_loci
#' @keywords internal
as_locus_windows <- function(windows, flank = 0) {
  if(is.character(windows)) {
    parsed <- lapply(windows, parse_locus, flank = flank)
    window_df <- data.frame(seqnames = vapply(parsed, `[[`, character(1), "chrom"),
                            start = vapply(parsed, `[[`, numeric(1), "start"),
                            end = vapply(parsed, `[[`, numeric(1), "end"))
    window_df$label <- if(!is.null(names(windows))) names(windows) else NA_character_
  } else {
    window_df <- as_genomic_df(windows)
    if(!is.data.frame(window_df) || !all(c("seqnames", "start", "end") %in% names(window_df))) {
      cli::cli_abort("{.arg windows} must be window strings, a GRanges, or a data frame with {.field seqnames}, {.field start}, and {.field end} columns.")
    }
    window_df <- data.frame(seqnames = as.character(window_df$seqnames),
                            start = window_df$start - flank, end = window_df$end + flank,
                            label = if("label" %in% names(window_df)) as.character(window_df$label) else NA_character_)
  }

  # Unlabeled windows are named by chromosome
  unlabeled <- is.na(window_df$label) | window_df$label == ""
  window_df$label[unlabeled] <- window_df$seqnames[unlabeled]
  window_df
}

#' @name locus_layout
#' @title Lay out genomic windows on a linear axis
#'
#' @description
#' Assigns each window a span of the linear axis. The first window starts at its own
#' genomic start, so a single window keeps genomic coordinates unchanged. For internal
#' use only.
#'
#' @param windows Data frame of windows from [as_locus_windows()]
#' @inheritParams coord_loci
#' @returns Data frame of windows with `chrom`, `chrom_key` (no `chr` prefix), `start`,
#'   `end`, `label`, `scale` (axis units per bp), `lin_start`, and `lin_end`
#' @keywords internal
locus_layout <- function(windows, widths = "proportional", gap = 0.02) {
  window_bp <- windows$end - windows$start
  axis_span <- if(widths == "equal") rep(mean(window_bp), length(window_bp)) else window_bp
  gap_span <- gap * sum(axis_span)

  lin_start <- windows$start[1] + c(0, cumsum(axis_span + gap_span)[-length(axis_span)])
  data.frame(chrom = windows$seqnames,
             chrom_key = stringr::str_remove(windows$seqnames, "^chr"),
             start = windows$start, end = windows$end, label = windows$label,
             scale = axis_span / window_bp,
             lin_start = lin_start, lin_end = lin_start + axis_span)
}

#' @name locus_breaks
#' @title Axis breaks within each window
#'
#' @description
#' Pretty breaks inside each window, in proportion to its share of the axis, labeled in
#' Mbp for windows wider than 10 kbp and bp otherwise. With several windows, breaks near
#' window edges are dropped so labels of neighboring windows do not collide. For
#' internal use only.
#'
#' @param layout Window layout from [locus_layout()]
#' @param n Approximate number of breaks across the whole axis, default: 8, or 5 for
#'   a single window
#' @returns List of break positions on the linear axis (`lin`) and their labels
#' @keywords internal
locus_breaks <- function(layout, n = 8) {
  share <- (layout$lin_end - layout$lin_start) / sum(layout$lin_end - layout$lin_start)
  # A single window keeps ggplot2's usual number of breaks and breaks at its edges
  if(nrow(layout) == 1) {
    n <- 5
  }
  edge_room <- if(nrow(layout) > 1) 0.1 else 0
  per_window <- lapply(seq_len(nrow(layout)), function(i) {
    window_breaks <- scales::breaks_pretty(n = max(2, round(n * share[i])))(c(layout$start[i], layout$end[i]))
    window_room <- edge_room * (layout$end[i] - layout$start[i])
    window_breaks <- window_breaks[window_breaks >= layout$start[i] + window_room &
                                     window_breaks <= layout$end[i] - window_room]
    if(length(window_breaks) == 0) {
      window_breaks <- round((layout$start[i] + layout$end[i]) / 2)
    }
    break_labels <- if(layout$end[i] - layout$start[i] + 1 > 10000) {
      scales::label_comma(accuracy = 0.001, scale = 1e-6)(window_breaks)
    } else {
      scales::label_comma(accuracy = 1)(window_breaks)
    }
    list(lin = to_linear(window_breaks, layout, i), label = break_labels)
  })
  list(lin = unlist(lapply(per_window, `[[`, "lin")),
       label = unlist(lapply(per_window, `[[`, "label")))
}

#' @name to_linear
#' @title Convert genomic positions in a window to the linear axis
#'
#' @description
#' For internal use only.
#'
#' @param position Genomic positions
#' @param layout Window layout from [locus_layout()]
#' @param window Row of the layout the positions belong to, one or one per position
#' @keywords internal
to_linear <- function(position, layout, window) {
  layout$lin_start[window] + (position - layout$start[window]) * layout$scale[window]
}

#' @name linearize_layer
#' @title Convert a layer's data to the linear axis
#'
#' @description
#' Converts `start`/`end` positions, or both anchors of two-anchor data, to the linear
#' axis of the window layout, and drops rows out of view. With more than one window,
#' ranges are clipped to window edges and split across windows. Data without position
#' columns but with `seqnames` is subset to the window chromosomes. For internal use
#' only.
#'
#' @param data Data frame of a single layer
#' @param layout Window layout from [locus_layout()]
#' @keywords internal
linearize_layer <- function(data, layout) {
  if(!is.data.frame(data)) {
    return(data)
  }
  columns <- names(data)

  if(all(c("seqnames1", "start1", "seqnames2", "start2") %in% columns)) {
    return(linearize_anchor_pairs(as.data.frame(data), layout))
  }
  if(!"seqnames" %in% columns) {
    return(data)
  }

  data <- as.data.frame(data)
  chrom_key <- chrom_keys(data$seqnames)
  if(!"start" %in% columns) {
    return(data[chrom_key %in% layout$chrom_key, , drop = FALSE])
  }

  # Rows overlapping each window, clipped to it when windows sit side by side
  has_end <- "end" %in% columns
  range_end <- if(has_end) data$end else data$start
  clip_ranges <- nrow(layout) > 1
  rows_by_chrom <- split(seq_len(nrow(data)), chrom_key)
  window_rows <- lapply(seq_len(nrow(layout)), function(i) {
    on_chrom <- rows_by_chrom[[layout$chrom_key[i]]]
    on_chrom[range_end[on_chrom] >= layout$start[i] & data$start[on_chrom] <= layout$end[i]]
  })

  # One subset for all windows, then positions converted with each row's window
  row_window <- rep(seq_len(nrow(layout)), lengths(window_rows))
  rows <- unlist(window_rows)
  linear_data <- data[rows, , drop = FALSE]
  row_start <- linear_data$start
  row_end <- range_end[rows]
  if(clip_ranges) {
    row_start <- pmax(row_start, layout$start[row_window])
    row_end <- pmin(row_end, layout$end[row_window])
  }
  linear_data$start <- to_linear(row_start, layout, row_window)
  if(has_end) {
    linear_data$end <- to_linear(row_end, layout, row_window)
  }
  linear_data
}

#' @name linearize_anchor_pairs
#' @title Convert both anchors of two-anchor data to the linear axis
#'
#' @description
#' Each anchor is converted using the first window it overlaps, or set to `NA` when it
#' is out of view. Rows with neither anchor in view are dropped. For internal use only.
#'
#' @inheritParams linearize_layer
#' @keywords internal
linearize_anchor_pairs <- function(data, layout) {
  in_view <- rep(FALSE, nrow(data))
  for(anchor in c("1", "2")) {
    anchor_start <- data[[paste0("start", anchor)]]
    end_column <- paste0("end", anchor)
    anchor_end <- if(end_column %in% names(data)) data[[end_column]] else anchor_start
    chrom_key <- chrom_keys(data[[paste0("seqnames", anchor)]])

    anchor_window <- rep(NA_integer_, nrow(data))
    for(i in seq_len(nrow(layout))) {
      in_window <- is.na(anchor_window) & chrom_key == layout$chrom_key[i] &
        anchor_end >= layout$start[i] & anchor_start <= layout$end[i]
      anchor_window[in_window] <- i
    }

    linear_start <- rep(NA_real_, nrow(data))
    linear_end <- rep(NA_real_, nrow(data))
    found <- which(!is.na(anchor_window))
    linear_start[found] <- to_linear(anchor_start[found], layout, anchor_window[found])
    linear_end[found] <- to_linear(anchor_end[found], layout, anchor_window[found])
    data[[paste0("start", anchor)]] <- linear_start
    if(end_column %in% names(data)) {
      data[[end_column]] <- linear_end
    }
    in_view <- in_view | !is.na(anchor_window)
  }
  data[in_view, , drop = FALSE]
}

#' @name chrom_keys
#' @title Chromosome names without the chr prefix
#'
#' @description
#' Strips the prefix from the unique names only, so it stays fast for millions of
#' rows. For internal use only.
#'
#' @param seqnames Chromosome names, as a character vector or factor
#' @keywords internal
chrom_keys <- function(seqnames) {
  if(is.factor(seqnames)) {
    return(stringr::str_remove(levels(seqnames), "^chr")[as.integer(seqnames)])
  }
  seqnames <- as.character(seqnames)
  unique_names <- unique(seqnames)
  stringr::str_remove(unique_names, "^chr")[match(seqnames, unique_names)]
}

#' @name parse_locus
#' @title Parse a genomic window string
#'
#' @description
#' Parses a window string and extends it by the flank on both ends. For internal use only.
#'
#' @inheritParams coord_locus
#' @returns List with `chrom`, `start`, and `end` of the locus
#' @keywords internal
parse_locus <- function(window_of_interest, flank = 0) {
  locus_gr <- gUtils::parse.gr(window_of_interest) + flank
  list("chrom" = as.character(GenomeInfoDb::seqnames(locus_gr)),
       "start" = GenomicRanges::start(locus_gr),
       "end" = GenomicRanges::end(locus_gr))
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
  on_chrom <- chrom_keys(data$seqnames) == stringr::str_remove(locus$chrom, "^chr")

  # Match overlap with window, using start alone for single position data
  in_window <- TRUE
  if("start" %in% names(data)) {
    range_end <- if("end" %in% names(data)) data$end else data$start
    in_window <- range_end >= locus$start & data$start <= locus$end
  }

  data[on_chrom & in_window, , drop = FALSE]
}
