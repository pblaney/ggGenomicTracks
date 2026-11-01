#' @name geom_ideogram
#' @title Chromosome ideogram
#'
#' @description
#' Draws the cytobands of a single chromosome between y = 0 and y = 1, colored by
#' Giemsa stain with a pinched centromere and an outline. Pair with
#' [scale_fill_giemsa()] for the standard band colors, [annotate_ideogram_window()]
#' and [annotate_ideogram_zoom()] to mark a locus, and [theme_ideogram()]. See
#' [plot_ideogram()] for all of these in one call.
#'
#' @param chrom Chromosome to draw, e.g. `"chr3"`.
#' @param data GRanges object or data frame of cytobands with `seqnames`, `start`,
#'   `end`, `name`, and `gieStain` columns. Defaults to [cytobands_hg38].
#' @param outline_colour Color of the chromosome outline.
#' @param linewidth Width of the chromosome outline.
#'
#' @returns List of ggplot2 layers
#' @export
#' @keywords geom
#'
#' @examples
#' \dontrun{
#' ggplot2::ggplot() +
#'   geom_ideogram("chr3") +
#'   annotate_ideogram_window("chr3:187734001-187751500") +
#'   scale_fill_giemsa() +
#'   theme_ideogram()
#' }
geom_ideogram <- function(chrom, data = NULL, outline_colour = "gray20", linewidth = 0.3) {
  bands <- get_ideogram_bands(chrom = chrom, data = data)
  chrom_start <- min(bands$start)
  chrom_end <- max(bands$end)

  # Centromere is drawn as two triangles meeting at the p/q boundary
  acen <- bands[bands$gieStain == "acen", ]
  if(nrow(acen) == 2) {
    cen_start <- min(acen$start)
    cen_mid <- acen$end[which.min(acen$start)]
    cen_end <- max(acen$end)
    centromere <- data.frame("x" = c(cen_start, cen_start, cen_mid, cen_mid, cen_end, cen_end),
                             "y" = c(0, 1, 0.5, 0.5, 1, 0),
                             "arm" = rep(c("p", "q"), each = 3))
    outline <- data.frame("x" = c(chrom_start, cen_start, cen_mid, cen_end, chrom_end,
                                  chrom_end, cen_end, cen_mid, cen_start, chrom_start),
                          "y" = c(0, 0, 0.5, 0, 0, 1, 1, 0.5, 1, 1))
  } else {
    # No centromere annotation, draw a plain bar
    centromere <- data.frame("x" = numeric(0), "y" = numeric(0), "arm" = character(0))
    outline <- data.frame("x" = c(chrom_start, chrom_end, chrom_end, chrom_start),
                          "y" = c(0, 0, 1, 1))
  }

  list(ggplot2::geom_rect(data = bands[bands$gieStain != "acen", ],
                          ggplot2::aes(xmin = start, xmax = end, ymin = 0, ymax = 1, fill = gieStain),
                          colour = NA, inherit.aes = FALSE),
       ggplot2::geom_polygon(data = centromere,
                             ggplot2::aes(x = x, y = y, group = arm),
                             fill = giemsa_colors[["acen"]], inherit.aes = FALSE),
       ggplot2::geom_polygon(data = outline,
                             ggplot2::aes(x = x, y = y),
                             fill = NA, colour = outline_colour, linewidth = linewidth, inherit.aes = FALSE))
}

#' @name scale_fill_giemsa
#' @title Giemsa stain fill scale for ideograms
#'
#' @description
#' Grayscale band colors by Giemsa stain intensity, with highlighted centromeres,
#' variable regions, and stalks.
#'
#' @param ... Other arguments passed to [ggplot2::scale_fill_manual()], e.g. `values`
#'   to override individual stain colors.
#'
#' @returns ggplot2 scale
#' @export
#' @keywords scale
scale_fill_giemsa <- function(...) {
  args <- utils::modifyList(list(values = giemsa_colors, guide = "none"), list(...))
  do.call(ggplot2::scale_fill_manual, args)
}

# Giemsa stain colors, from unstained (gneg) to darkest (gpos100)
giemsa_colors <- c("gneg" = "#FFFFFF", "gpos25" = "#C8C8C8", "gpos50" = "#969696",
                   "gpos75" = "#5A5A5A", "gpos100" = "#1E1E1E", "gvar" = "#DCDCDC",
                   "stalk" = "#8296B4", "acen" = "#C0504D")

#' @name annotate_ideogram_window
#' @title Highlight a window on an ideogram
#'
#' @description
#' Draws a box over the window on a [geom_ideogram()]. The outline keeps the box visible
#' even when the window is far narrower than a pixel at chromosome scale.
#'
#' @inheritParams coord_locus
#' @param fill Fill color of the box.
#' @param colour Outline color of the box.
#' @param linewidth Outline width of the box.
#'
#' @returns ggplot2 layer
#' @export
#' @keywords annotation
annotate_ideogram_window <- function(window_of_interest, flank = 0, fill = "#FFA319",
                                     colour = fill, linewidth = 1) {
  locus <- parse_locus(window_of_interest = window_of_interest, flank = flank)
  ggplot2::annotate("rect", xmin = locus$start, xmax = locus$end, ymin = -0.1, ymax = 1.1,
                    fill = fill, colour = colour, linewidth = linewidth)
}

#' @name annotate_ideogram_zoom
#' @title Zoom trace from an ideogram window to the track below
#'
#' @description
#' Draws a funnel from the window on a [geom_ideogram()] down to the full chromosome
#' width at the bottom of the panel. When the ideogram is stacked above locus tracks
#' with patchwork, the funnel connects the window to the zoomed tracks.
#'
#' @inheritParams coord_locus
#' @inheritParams geom_ideogram
#' @param fill Fill color of the funnel.
#' @param alpha Transparency of the funnel.
#'
#' @returns ggplot2 layer
#' @export
#' @keywords annotation
annotate_ideogram_zoom <- function(window_of_interest, flank = 0, data = NULL,
                                   fill = "#FFA319", alpha = 0.35) {
  locus <- parse_locus(window_of_interest = window_of_interest, flank = flank)
  bands <- get_ideogram_bands(chrom = locus$chrom, data = data)
  ggplot2::annotate("polygon",
                    x = c(locus$start, locus$end, max(bands$end), min(bands$start)),
                    y = c(-0.1, -0.1, -1.2, -1.2),
                    fill = fill, alpha = alpha)
}

#' @name theme_ideogram
#' @title Theme for ideogram plots
#'
#' @description
#' Blank theme with a centered title and chromosome position labels on the x-axis.
#'
#' @param base_size Base font size, default: 11
#'
#' @returns ggplot2 theme object
#' @export
#' @keywords theme
theme_ideogram <- function(base_size = 11) {
  ggplot2::theme_void(base_size = base_size) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = ggplot2::rel(1.1), hjust = 0.5),
                   axis.text.x = ggplot2::element_text(size = ggplot2::rel(0.82)),
                   axis.ticks.x = ggplot2::element_line(color = "black"),
                   axis.ticks.length.x = ggplot2::unit(1.5, units = "mm"),
                   plot.margin = ggplot2::unit(c(2, 2, 0.5, 2), "mm"))
}

#' @name plot_ideogram
#' @title Ideogram of a locus' chromosome with zoom trace
#'
#' @description
#' Complete ideogram plot for the top of a track stack: the chromosome with the window
#' highlighted, a zoom trace down to the tracks below, chromosome positions on top,
#' and a title of the chromosome and cytoband span, e.g. "Chr 3 (q27.3)". Built from
#' [geom_ideogram()] and its annotation layers, so it can be modified with `+`.
#'
#' @inheritParams coord_locus
#' @inheritParams geom_ideogram
#' @param highlight_fill Color of the window highlight and zoom trace.
#' @param zoom_trace Draw the zoom trace, default: TRUE
#'
#' @returns ggplot object
#' @export
#' @keywords plot
#'
#' @examples
#' \dontrun{
#' win <- "chr3:187734001-187751500"
#' patchwork::wrap_plots(plot_ideogram(win, flank = 5000),
#'                       bed_plot + coord_locus(win, flank = 5000),
#'                       ncol = 1, heights = c(0.4, 1))
#' }
plot_ideogram <- function(window_of_interest, flank = 0, data = NULL,
                          highlight_fill = "#FFA319", zoom_trace = TRUE) {
  locus <- parse_locus(window_of_interest = window_of_interest, flank = flank)

  ideogram <- ggplot2::ggplot()
  if(zoom_trace) {
    ideogram <- ideogram +
      annotate_ideogram_zoom(window_of_interest = window_of_interest, flank = flank,
                             data = data, fill = highlight_fill)
  }

  ideogram +
    geom_ideogram(chrom = locus$chrom, data = data) +
    annotate_ideogram_window(window_of_interest = window_of_interest, flank = flank, fill = highlight_fill) +
    scale_fill_giemsa() +
    ggplot2::scale_x_continuous(position = "top",
                                labels = scales::label_comma(scale = 1e-6, suffix = " Mb"),
                                expand = c(0.01, 0.01)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    # An explicit coord marks the ideogram as whole-chromosome, so stack_tracks() leaves it
    ggplot2::coord_cartesian() +
    ggplot2::labs(title = locus_cytoband_label(window_of_interest = window_of_interest,
                                               flank = flank, data = data)) +
    theme_ideogram()
}

#' @name get_ideogram_bands
#' @title Cytobands of a single chromosome
#'
#' @description
#' Subsets cytoband data to one chromosome, sorted by position. For internal use only.
#'
#' @inheritParams geom_ideogram
#' @keywords internal
get_ideogram_bands <- function(chrom, data = NULL) {
  if(is.null(data)) {
    data <- get_data(name_of_data = "cytobands_hg38")
  }
  bands <- as_genomic_df(data)
  bands <- bands[as.character(bands$seqnames) == chrom, ]
  if(nrow(bands) == 0) {
    cli::cli_abort("No cytobands found for chromosome {.val {chrom}}.")
  }
  bands[order(bands$start), ]
}

#' @name locus_cytoband_label
#' @title Chromosome and cytoband span label for a locus
#'
#' @description
#' Builds a label like "Chr 3 (q27.3)" or "Chr 3 (q27.2-q27.3)" from the cytobands
#' overlapping the window. For internal use only.
#'
#' @inheritParams coord_locus
#' @inheritParams geom_ideogram
#' @keywords internal
locus_cytoband_label <- function(window_of_interest, flank = 0, data = NULL) {
  locus <- parse_locus(window_of_interest = window_of_interest, flank = flank)
  overlap <- subset_to_locus(get_ideogram_bands(chrom = locus$chrom, data = data), locus = locus)
  band_span <- unique(c(overlap$name[1], overlap$name[nrow(overlap)]))
  paste0(stringr::str_replace(string = locus$chrom, pattern = "^chr", replacement = "Position on chromosome "),
         " (", paste(band_span, collapse = "-"), ")")
}
