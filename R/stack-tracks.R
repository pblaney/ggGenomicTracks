#' @name stack_tracks
#' @title Stack tracks over shared genomic windows
#'
#' @description
#' Stacks tracks vertically with patchwork, giving every track the same genomic
#' windows so their x-axes line up. Each track gets the same [coord_loci()] or
#' [coord_genome()], window labels are kept on the top track only, and the x-axis
#' labels on the bottom track only. Tracks with another coord, such as the ideogram
#' from [plot_ideogram()], are stacked as they are.
#'
#' @param ... Tracks to stack, top to bottom, as ggplots or a single list of ggplots.
#'   Tracks need no coord of their own.
#' @param windows Windows to show on every track, as for [coord_loci()], or
#'   `"genome"` for [coord_genome()]. Default: the windows of the first track that
#'   has a [coord_loci()] or [coord_locus()].
#' @param flank Number of bp to extend both ends of each window, default: 0
#' @param heights Relative heights of the tracks, default: equal
#' @param upper_x_axis Bottom x-axis style of the tracks above the bottom track:
#'   `"ticks"` (default), `"blank"`, or `"keep"` to leave each track's own theme.
#' @param guides Passed to [patchwork::wrap_plots()], default: `"collect"` to gather
#'   legends on the right.
#' @param highlight Regions to highlight across the stack, as for [highlight_stack()].
#'   Default: none.
#' @param highlight_args List of other arguments passed to [highlight_stack()], e.g.
#'   `list(fill = c("#9BD08A", "#DC766F"), legend_title = "Windows")`.
#' @inheritParams coord_loci
#'
#' @returns A patchwork, or with `highlight`, a highlighted stack from
#'   [highlight_stack()]
#' @export
#' @keywords annotation
#'
#' @examples
#' \dontrun{
#' windows <- c(MYC = "chr8:126500001-129500000", IGH = "chr14:105550001-106100000")
#'
#' stack_tracks(
#'   ggplot2::ggplot(junctions) + geom_junction(ggplot2::aes(colour = type)) + theme_genomic(),
#'   ggplot2::ggplot(cn) + geom_copy_number(cn_cap = 8) + theme_genomic(),
#'   ggplot2::ggplot() + geom_gene_track(genes_of_interest = c("MYC", "IGHM")) +
#'     theme_genomic(y_axis = "blank"),
#'   windows = windows, heights = c(1, 1, 0.5))
#' }
stack_tracks <- function(..., windows = NULL, flank = 0, widths = c("proportional", "equal"), gap = 0.02,
                         dividers = TRUE, window_labels = NULL, heights = NULL,
                         upper_x_axis = c("ticks", "blank", "keep"), guides = "collect",
                         highlight = NULL, highlight_args = list()) {
  widths <- match.arg(widths)
  upper_x_axis <- match.arg(upper_x_axis)
  tracks <- list(...)
  if(length(tracks) == 1 && is.list(tracks[[1]]) && !ggplot2::is_ggplot(tracks[[1]])) {
    tracks <- tracks[[1]]
  }

  # Locus tracks use coord_loci() or no coord of their own; others are left as they are
  is_locus_track <- vapply(tracks, function(track) {
    ggplot2::is_ggplot(track) && (inherits(track$coordinates, "CoordLoci") || isTRUE(track$coordinates$default))
  }, logical(1))
  locus_index <- which(is_locus_track)
  if(length(locus_index) == 0) {
    cli::cli_abort("None of the tracks can take shared windows; they all have a coord other than {.fn coord_loci}.")
  }

  # Shared coord, from the windows given or from the first track with windows
  if(identical(windows, "genome")) {
    coord_spec <- list(fun = "coord_genome", args = list())
  } else if(!is.null(windows)) {
    coord_spec <- list(fun = "coord_loci",
                       args = list(windows = windows, flank = flank, widths = widths, gap = gap, dividers = dividers))
  } else {
    with_windows <- Filter(function(track) inherits(track$coordinates, "CoordLoci"), tracks[locus_index])
    if(length(with_windows) == 0) {
      cli::cli_abort("Give {.arg windows}, or add {.fn coord_loci} to one of the tracks.")
    }
    coord_spec <- with_windows[[1]]$coordinates$spec
  }

  # Window labels on the top track and x-axis labels on the bottom track only
  for(i in locus_index) {
    coord_args <- coord_spec$args
    if(coord_spec$fun == "coord_loci") {
      coord_args$window_labels <- if(i == locus_index[1]) window_labels else FALSE
    }
    tracks[[i]] <- suppressMessages(tracks[[i]] + do.call(coord_spec$fun, coord_args))
    if(i != utils::tail(locus_index, 1) && upper_x_axis != "keep") {
      tracks[[i]] <- tracks[[i]] + upper_axis_theme(upper_x_axis)
    }
  }

  stack <- patchwork::wrap_plots(tracks, ncol = 1, heights = heights, guides = guides)
  if(is.null(highlight)) {
    return(stack)
  }
  do.call(highlight_stack, c(list(stack = stack, regions = highlight), highlight_args))
}

#' @name upper_axis_theme
#' @title Bottom x-axis theme of upper tracks in a stack
#'
#' @description
#' For internal use only.
#'
#' @param style `"ticks"` or `"blank"`
#' @keywords internal
upper_axis_theme <- function(style) {
  switch(style,
    "ticks" = ggplot2::theme(axis.text.x.bottom = ggplot2::element_blank(),
                             axis.line.x.bottom = ggplot2::element_blank(),
                             axis.ticks.x.bottom = ggplot2::element_line(colour = "black"),
                             axis.ticks.length.x.bottom = ggplot2::unit(1.25, "mm")),
    "blank" = ggplot2::theme(axis.text.x.bottom = ggplot2::element_blank(),
                             axis.line.x.bottom = ggplot2::element_blank(),
                             axis.ticks.x.bottom = ggplot2::element_blank()))
}
