#' @name highlight_stack
#' @title Highlight sub-regions across a stack of tracks
#'
#' @description
#' Shades one or more sub-regions as unbroken bands running through every track of a
#' stack, including the space between tracks. Bands cover each run of adjacent tracks
#' that share the same [coord_locus()] window; other plots in the stack, such as the
#' ideogram from [plot_ideogram()], are left alone. Apply it as the last step, after
#' the stack is assembled with patchwork.
#'
#' Bands are drawn over the data by default, so they stay visible on tracks whose data
#' covers the whole panel, like dense bed tracks. Adjust `alpha` to set how strongly
#' they tint the data, or use `layer = "below"` to keep the data untinted.
#'
#' @param stack A patchwork of tracks, e.g. from `patchwork::wrap_plots(..., ncol = 1)`,
#'   or a single ggplot.
#' @param regions Regions to highlight as one of:
#'   * A character vector of windows, e.g. `"chr3:187744872-187744917"`. Names, if
#'     given, are used as labels.
#'   * A data frame with `start` and `end` columns and optional `seqnames`, `label`,
#'     and `fill` columns.
#'   * A GRanges object, with optional `label` and `fill` metadata columns.
#' @param fill Fill colors, either one color for all regions or one per label (named
#'   or in order of appearance). Ignored if `regions` has a `fill` column. By default a
#'   single region is orange and labeled regions use a hue palette.
#' @param alpha Opacity of the bands, default: 0.3
#' @param colour Outline color of the bands, default: none
#' @param layer Draw the bands `"above"` (default) or `"below"` the track data.
#' @param tracks Positions of the plots in the stack to cover, e.g. `3:5`. By default
#'   all plots using [coord_locus()]. Bands still run unbroken through adjacent tracks.
#' @param show_legend Add a legend of the region labels, default: TRUE when regions are
#'   labeled.
#' @param legend_title Title of the legend, default: "Highlight"
#'
#' @returns A gtable of the highlighted stack with class `highlighted_stack`. Print it to
#'   draw, save it with [ggplot2::ggsave()], or wrap it with
#'   `patchwork::wrap_elements(full = )` to combine with other plots.
#' @export
#' @keywords annotation
#'
#' @examples
#' \dontrun{
#' win <- "chr3:187744850-187745350"
#'
#' stack <- patchwork::wrap_plots(
#'   plot_ideogram(win, flank = 250),
#'   ggplot2::ggplot(snvs) + geom_count_mutations(...) + coord_locus(win, 250) +
#'     theme_genomic(x_axis = "ticks"),
#'   ggplot2::ggplot() + geom_gene_track(genes_of_interest = "BCL6") +
#'     coord_locus(win, 250) + theme_genomic(y_axis = "blank"),
#'   ncol = 1, guides = "collect")
#'
#' highlight_stack(stack,
#'                 c("Normal MBC Hotspot 1" = "chr3:187744872-187744917",
#'                   "Normal MBC Hotspot 2" = "chr3:187745199-187745246",
#'                   "DLBCL Hotspot" = "chr3:187744678-187744698"),
#'                 fill = c("#9BD08A", "#C2EFB4", "#DC766F"),
#'                 legend_title = "Windows")
#' }
highlight_stack <- function(stack, regions, fill = NULL, alpha = 0.3, colour = NA,
                            layer = c("above", "below"), tracks = NULL, show_legend = NULL,
                            legend_title = "Highlight") {
  layer <- match.arg(layer)
  regions <- as_highlight_regions(regions, fill = fill)
  if(is.null(show_legend)) {
    show_legend <- any(!is.na(regions$label))
  }

  # Lay out the stack and find each plot's panel
  plots <- stack_plots(stack)
  if(inherits(stack, "patchwork")) {
    stack_gtable <- patchwork::patchworkGrob(stack)
  } else {
    stack_gtable <- ggplot2::ggplotGrob(stack)
  }
  panel_layout <- stack_gtable$layout[grepl("^panel(-[0-9]+)?$", stack_gtable$layout$name), , drop = FALSE]
  if(nrow(panel_layout) != length(plots)) {
    cli::cli_abort(c("Each plot in {.arg stack} must have exactly one panel.",
                     "i" = "Found {nrow(panel_layout)} panel{?s} for {length(plots)} plot{?s}; faceted or nested plots are not supported."))
  }
  panel_order <- as.integer(stringr::str_extract(panel_layout$name, "[0-9]+$"))
  panel_layout <- panel_layout[order(panel_order), , drop = FALSE]

  # Window of each locus track, as the chromosome and visible x-range, NA otherwise
  track_windows <- vapply(plots, function(plot) {
    if(!inherits(plot$coordinates, "CoordLocus")) {
      return(NA_character_)
    }
    x_range <- ggplot2::ggplot_build(plot)$layout$panel_params[[1]]$x.range
    paste(plot$coordinates$locus$chrom, x_range[1], x_range[2])
  }, character(1))
  if(!is.null(tracks)) {
    not_locus <- tracks[is.na(track_windows[tracks])]
    if(length(not_locus) > 0) {
      cli::cli_warn("Skipping plot{?s} {not_locus} in {.arg tracks}, which do{?es/} not use {.fn coord_locus}.")
    }
    track_windows[setdiff(seq_along(track_windows), tracks)] <- NA
  }
  if(all(is.na(track_windows))) {
    cli::cli_abort("No plots to highlight; none of the selected plots in {.arg stack} use {.fn coord_locus}.")
  }

  # One band per region for each run of adjacent tracks sharing a window
  window_runs <- rle(track_windows)
  run_ends <- cumsum(window_runs$lengths)
  run_starts <- run_ends - window_runs$lengths + 1
  for(run in which(!is.na(window_runs$values))) {
    window_parts <- strsplit(window_runs$values[run], " ")[[1]]
    run_regions <- regions[is.na(regions$seqnames) |
                             stringr::str_remove(regions$seqnames, "^chr") == stringr::str_remove(window_parts[1], "^chr"), ,
                           drop = FALSE]
    if(nrow(run_regions) == 0) {
      next
    }
    x_range <- as.numeric(window_parts[2:3])
    run_panels <- panel_layout[run_starts[run]:run_ends[run], , drop = FALSE]

    band_grob <- grid::rectGrob(x = grid::unit((run_regions$start - x_range[1]) / diff(x_range), "npc"),
                                width = grid::unit((run_regions$end - run_regions$start) / diff(x_range), "npc"),
                                y = grid::unit(0, "npc"), height = grid::unit(1, "npc"),
                                just = c("left", "bottom"),
                                gp = grid::gpar(fill = scales::alpha(run_regions$fill, alpha), col = colour))
    stack_gtable <- gtable::gtable_add_grob(stack_gtable, band_grob,
                                            t = min(run_panels$t), b = max(run_panels$b),
                                            l = min(run_panels$l), r = max(run_panels$r),
                                            z = if(layer == "below") 0.5 else Inf, clip = "on",
                                            name = paste0("highlight-", run))
  }

  # Legend of labeled regions in a column to the right of the stack
  if(show_legend && any(!is.na(regions$label))) {
    stack_gtable <- add_highlight_legend(stack_gtable, regions, alpha = alpha, colour = colour,
                                         legend_title = legend_title)
  }

  class(stack_gtable) <- c("highlighted_stack", class(stack_gtable))
  stack_gtable
}

#' @export
print.highlighted_stack <- function(x, newpage = TRUE, ...) {
  if(newpage) {
    grid::grid.newpage()
  }
  grid::grid.draw(x)
  invisible(x)
}

#' @export
plot.highlighted_stack <- print.highlighted_stack

#' @name as_highlight_regions
#' @title Standardize regions to highlight
#'
#' @description
#' Converts window strings, data frames, or GRanges to a data frame with `seqnames`,
#' `start`, `end`, `label`, and `fill` columns. For internal use only.
#'
#' @inheritParams highlight_stack
#' @keywords internal
as_highlight_regions <- function(regions, fill = NULL) {
  if(is.character(regions)) {
    windows <- lapply(regions, parse_locus)
    regions <- data.frame(seqnames = vapply(windows, `[[`, character(1), "chrom"),
                          start = vapply(windows, `[[`, numeric(1), "start"),
                          end = vapply(windows, `[[`, numeric(1), "end"))
    if(!is.null(names(windows))) {
      regions$label <- names(windows)
    }
  }
  regions <- as_genomic_df(regions)
  if(!is.data.frame(regions) || !all(c("start", "end") %in% names(regions))) {
    cli::cli_abort("{.arg regions} must be window strings, a GRanges, or a data frame with {.field start} and {.field end} columns.")
  }

  regions$seqnames <- if("seqnames" %in% names(regions)) as.character(regions$seqnames) else NA_character_
  regions$label <- if("label" %in% names(regions)) as.character(regions$label) else NA_character_
  regions$label[regions$label == ""] <- NA

  # Fill per region from a fill column, one color, or one color per label
  if(!"fill" %in% names(regions)) {
    labels <- unique(regions$label[!is.na(regions$label)])
    if(is.null(fill)) {
      fill <- if(length(labels) <= 1) "#FFA319" else scales::pal_hue()(length(labels))
    }
    if(length(fill) == 1) {
      regions$fill <- fill
    } else {
      if(is.null(names(fill))) {
        if(length(fill) != length(labels)) {
          cli::cli_abort("{.arg fill} must be one color or one per label ({length(labels)}), not {length(fill)}.")
        }
        names(fill) <- labels
      }
      regions$fill <- unname(fill[regions$label])
    }
  }
  regions[, c("seqnames", "start", "end", "label", "fill")]
}

#' @name stack_plots
#' @title List the plots of a patchwork in layout order
#'
#' @description
#' For internal use only.
#'
#' @inheritParams highlight_stack
#' @keywords internal
stack_plots <- function(stack) {
  if(!inherits(stack, "patchwork")) {
    return(list(stack))
  }
  # The patchwork object is its own last plot, with the others stored in patches
  last_plot <- stack
  last_plot$patches <- NULL
  class(last_plot) <- setdiff(class(last_plot), "patchwork")
  plots <- c(stack$patches$plots, list(last_plot))
  if(any(vapply(plots, inherits, logical(1), what = "patchwork"))) {
    cli::cli_abort("Nested patchworks are not supported in {.arg stack}; build it with a single {.fn patchwork::wrap_plots} call.")
  }
  plots
}

#' @name add_highlight_legend
#' @title Add a legend of highlighted regions to a stack
#'
#' @description
#' Builds the legend from a throwaway plot of the regions and adds it below the
#' collected legends of the stack, or in a new column to the right if legends are not
#' collected. For internal use only.
#'
#' @param stack_gtable Gtable of the stack
#' @param regions Standardized regions from [as_highlight_regions()]
#' @inheritParams highlight_stack
#' @keywords internal
add_highlight_legend <- function(stack_gtable, regions, alpha, colour, legend_title) {
  labeled <- regions[!is.na(regions$label), , drop = FALSE]
  labeled <- labeled[!duplicated(labeled$label), , drop = FALSE]
  labeled$label <- factor(labeled$label, levels = labeled$label)

  legend_plot <- ggplot2::ggplot(labeled) +
    ggplot2::geom_rect(ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = label),
                       alpha = alpha, colour = colour) +
    ggplot2::scale_fill_manual(name = legend_title, values = stats::setNames(labeled$fill, labeled$label)) +
    theme_genomic()
  legend_gtable <- ggplot2::ggplotGrob(legend_plot)
  legend_box <- legend_gtable$grobs[[which(grepl("^guide-box", legend_gtable$layout$name) &
                                             !vapply(legend_gtable$grobs, inherits, logical(1), what = "zeroGrob"))[1]]]

  # Without collected guides, give the legend its own column
  box_index <- which(stack_gtable$layout$name == "guide-box")
  if(length(box_index) == 0 || inherits(stack_gtable$grobs[[box_index]], "zeroGrob")) {
    stack_gtable <- gtable::gtable_add_cols(stack_gtable, gtable::gtable_width(legend_box))
    return(gtable::gtable_add_grob(stack_gtable, legend_box, t = 1, b = nrow(stack_gtable),
                                   l = ncol(stack_gtable), clip = "off", name = "highlight-legend"))
  }

  # Otherwise add it as a new row below the collected legends, left aligned
  guide_box <- stack_gtable$grobs[[box_index]]
  legend_rows <- guide_box$layout[guide_box$layout$name == "guides", , drop = FALSE]
  last_row <- max(legend_rows$b)
  legend_column <- legend_rows$l[1]
  highlight <- legend_box$grobs[[which(legend_box$layout$name == "guides")[1]]]
  box_width <- max(guide_box$widths[legend_column], gtable::gtable_width(highlight))
  aligned <- gtable::gtable(widths = grid::unit.c(gtable::gtable_width(highlight),
                                                  box_width - gtable::gtable_width(highlight)),
                            heights = gtable::gtable_height(highlight))
  aligned <- gtable::gtable_add_grob(aligned, highlight, t = 1, l = 1, clip = "off")

  guide_box <- gtable::gtable_add_rows(guide_box, grid::unit(11, "pt"), pos = last_row)
  guide_box <- gtable::gtable_add_rows(guide_box, gtable::gtable_height(highlight), pos = last_row + 1)
  guide_box <- gtable::gtable_add_grob(guide_box, aligned, t = last_row + 2, l = legend_column,
                                       clip = "off", name = "guides")
  guide_box$widths[legend_column] <- box_width
  stack_gtable$grobs[[box_index]] <- guide_box

  # Widen the guide column if the highlight legend is wider than the collected ones
  box_column <- stack_gtable$layout$l[box_index]
  stack_gtable$widths[box_column] <- max(stack_gtable$widths[box_column], gtable::gtable_width(guide_box))
  stack_gtable
}
