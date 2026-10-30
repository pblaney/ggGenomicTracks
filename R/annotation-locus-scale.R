#' @name annotation_locus_scale
#' @title Locus coordinate and scale bar annotation
#'
#' @description
#' Adds a scale bar and the window coordinates to a track panel. The bar length is
#' chosen from the visible x-range (a round 1, 2, or 5 x 10^n bp up to 20% of the
#' width), so it adapts to any zoom. With [coord_locus()] the coordinates show the
#' exact window, e.g. "chr3:187,729,001-187,756,500 (27,500 bp)"; otherwise the
#' visible range is shown without a chromosome.
#'
#' Add to any track, or recreate a standalone scale track with
#' `ggplot() + annotation_locus_scale(location = "middle") + coord_locus(win) + theme_void()`.
#'
#' @param location Vertical placement in the panel, one of `"top"`, `"middle"`, or `"bottom"`.
#' @param bar_side Side of the panel for the scale bar, `"right"` or `"left"`. The
#'   coordinates are placed on the opposite side.
#' @param show_coordinates Show the window coordinates, default: TRUE
#' @param colour Color of the bar and text.
#' @param linewidth Width of the bar in mm, default: 0.8
#' @param text_size Size of the text in mm, default: 3.4
#' @param pad Distance from the panel edges as a fraction of panel width, default: 0.02
#'
#' @returns ggplot2 layer
#' @export
#' @keywords annotation
#'
#' @examples
#' \dontrun{
#' ggplot2::ggplot() +
#'   geom_bed_track(data = chrom_states, ggplot2::aes(y = sample, fill = segment_label)) +
#'   annotation_locus_scale() +
#'   coord_locus("chr3:187734001-187751500", flank = 5000) +
#'   theme_genomic()
#' }
annotation_locus_scale <- function(location = c("top", "middle", "bottom"), bar_side = c("right", "left"),
                                   show_coordinates = TRUE, colour = "black", linewidth = 0.8,
                                   text_size = 3.4, pad = 0.02) {
  ggplot2::layer(data = data.frame(x = NA), mapping = NULL,
                 stat = ggplot2::StatIdentity, geom = GeomLocusScale, position = ggplot2::PositionIdentity,
                 inherit.aes = FALSE, show.legend = FALSE,
                 params = list(location = match.arg(location), bar_side = match.arg(bar_side),
                               show_coordinates = show_coordinates, colour = colour,
                               linewidth = linewidth, text_size = text_size, pad = pad))
}

#' @rdname annotation_locus_scale
#' @format NULL
#' @usage NULL
#' @export
GeomLocusScale <- ggplot2::ggproto("GeomLocusScale", ggplot2::Geom,
  draw_key = ggplot2::draw_key_blank,

  draw_panel = function(data, panel_params, coord, location = "top", bar_side = "right",
                        show_coordinates = TRUE, colour = "black", linewidth = 0.8,
                        text_size = 3.4, pad = 0.02) {
    # Visible x-range, and the exact window if zoomed with coord_locus()
    visible_range <- panel_params$x.range
    if(inherits(coord, "CoordLocus")) {
      window_start <- coord$locus$start
      window_end <- coord$locus$end
      window_chrom <- coord$locus$chrom
    } else {
      window_start <- ceiling(visible_range[1])
      window_end <- floor(visible_range[2])
      window_chrom <- NULL
    }

    # Bar spans a round length of up to 20% of the visible width
    bar_bp <- nice_scale_length(diff(visible_range) * 0.2)
    bar_npc <- bar_bp / diff(visible_range)
    if(bar_side == "right") {
      bar_x <- c(1 - pad - bar_npc, 1 - pad)
      coordinates_x <- pad
      coordinates_hjust <- 0
    } else {
      bar_x <- c(pad, pad + bar_npc)
      coordinates_x <- 1 - pad
      coordinates_hjust <- 1
    }

    # Bar with its label above and the two-line coordinates form one block, sized to
    # the taller of the two (in mm) and placed within the panel
    tick_half_mm <- 1.25
    label_gap_mm <- 0.75
    bar_group_mm <- 2 * tick_half_mm + label_gap_mm + text_size
    coordinates_mm <- if(show_coordinates) 2 * 1.2 * text_size else 0
    block_half <- grid::unit(max(bar_group_mm, coordinates_mm) / 2, "mm")
    edge_gap <- grid::unit(1.5, "mm")
    block_y <- switch(location,
                      "top" = grid::unit(1, "npc") - edge_gap - block_half,
                      "middle" = grid::unit(0.5, "npc"),
                      "bottom" = grid::unit(0, "npc") + edge_gap + block_half)
    bar_y <- block_y - grid::unit((label_gap_mm + text_size) / 2, "mm")
    tick_half <- grid::unit(tick_half_mm, "mm")

    line_gp <- grid::gpar(col = colour, lwd = linewidth * ggplot2::.pt, lineend = "butt")
    text_gp <- grid::gpar(col = colour, fontsize = text_size * ggplot2::.pt)

    scale_grobs <- grid::gList(
      grid::segmentsGrob(x0 = grid::unit(bar_x[1], "npc"), x1 = grid::unit(bar_x[2], "npc"),
                         y0 = bar_y, y1 = bar_y, gp = line_gp),
      grid::segmentsGrob(x0 = grid::unit(bar_x, "npc"), x1 = grid::unit(bar_x, "npc"),
                         y0 = bar_y - tick_half, y1 = bar_y + tick_half, gp = line_gp),
      grid::textGrob(label = format_bp(bar_bp), x = grid::unit(mean(bar_x), "npc"),
                     y = bar_y + tick_half + grid::unit(label_gap_mm, "mm"), vjust = 0, gp = text_gp)
    )

    if(show_coordinates) {
      coordinates_label <- paste0(if(!is.null(window_chrom)) paste0(window_chrom, ":"),
                                  scales::label_comma(accuracy = 1)(window_start), "-",
                                  scales::label_comma(accuracy = 1)(window_end),
                                  "\n(", scales::label_comma(accuracy = 1)(window_end - window_start + 1), " bp)")
      scale_grobs <- grid::gList(scale_grobs,
                                 grid::textGrob(label = coordinates_label, x = grid::unit(coordinates_x, "npc"),
                                                y = block_y, hjust = coordinates_hjust, gp = text_gp))
    }

    grid::gTree(children = scale_grobs, name = "locus_scale")
  }
)

#' @name nice_scale_length
#' @title Round scale bar length
#'
#' @description
#' Largest 1, 2, or 5 x 10^n value not exceeding the target length. For internal use only.
#'
#' @param target_bp Maximum scale bar length in bp
#' @keywords internal
nice_scale_length <- function(target_bp) {
  magnitude <- 10^floor(log10(target_bp))
  candidates <- magnitude * c(1, 2, 5)
  max(candidates[candidates <= target_bp])
}

#' @name format_bp
#' @title Format a length in bp, kb, or Mb
#'
#' @description
#' For internal use only.
#'
#' @param bp Length in bp
#' @keywords internal
format_bp <- function(bp) {
  if(bp >= 1e6) {
    paste0(bp / 1e6, " Mb")
  } else if(bp >= 1e3) {
    paste0(bp / 1e3, " kb")
  } else {
    paste0(bp, " bp")
  }
}
