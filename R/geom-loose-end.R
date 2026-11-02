#' @name geom_loose_end
#' @title Loose end track
#'
#' @description
#' Draws loose ends, breakpoints with copy number change but no junction partner as
#' inferred by JaBbA, as a short stem topped with an open circle. `x` defaults to the
#' `start` column. Map `y` to the copy number at each loose end, e.g. from
#' [loose_end_copy_number()], to place them on a copy number track; otherwise they sit
#' at the bottom of the panel.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()].
#' @param data Data frame of loose ends with `seqnames` and `start` columns, e.g. a
#'   JaBbA `loose_ends.bed`.
#' @param stem_height Height of the stems in mm, default: 3
#' @param ... Other arguments passed to [ggplot2::layer()], e.g. `colour` or
#'   `linewidth`.
#' @param na.rm,show.legend,inherit.aes See [ggplot2::layer()].
#'
#' @returns ggplot2 layer
#' @export
#' @keywords geom
#'
#' @examples
#' \dontrun{
#' loose <- data.table::fread("sample.jabba.phased.loose_ends.bed",
#'                            col.names = c("seqnames", "start", "end"))
#'
#' ggplot2::ggplot() +
#'   geom_copy_number(data = cn) +
#'   geom_loose_end(data = loose_end_copy_number(loose, cn), ggplot2::aes(y = cn)) +
#'   coord_genome() +
#'   theme_genomic()
#' }
geom_loose_end <- function(mapping = NULL, data = NULL, stem_height = 3, ...,
                           na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  loose_mapping <- ggplot2::aes(x = start, y = 0)
  for(aes_name in names(mapping)) {
    loose_mapping[[aes_name]] <- mapping[[aes_name]]
  }

  ggplot2::layer(geom = GeomLooseEnd, stat = ggplot2::StatIdentity, position = ggplot2::PositionIdentity,
                 data = as_genomic_df(data), mapping = loose_mapping,
                 show.legend = show.legend, inherit.aes = inherit.aes,
                 params = list(stem_height = stem_height, na.rm = na.rm, ...))
}

#' @rdname geom_loose_end
#' @format NULL
#' @usage NULL
#' @export
GeomLooseEnd <- ggplot2::ggproto("GeomLooseEnd", ggplot2::Geom,
  required_aes = c("x", "y"),
  default_aes = ggplot2::aes(colour = "#1F3A93", linewidth = 0.5, alpha = NA),
  draw_key = ggplot2::draw_key_path,

  draw_panel = function(data, panel_params, coord, stem_height = 3, na.rm = FALSE) {
    data <- coord$transform(data, panel_params)
    loose_gp <- grid::gpar(col = ggplot2::alpha(data$colour, data$alpha), lwd = data$linewidth * ggplot2::.pt,
                           fill = NA)
    stem_top <- grid::unit(data$y, "npc") + grid::unit(stem_height, "mm")
    grid::gTree(children = grid::gList(
      grid::segmentsGrob(x0 = data$x, x1 = data$x, y0 = grid::unit(data$y, "npc"), y1 = stem_top, gp = loose_gp),
      grid::circleGrob(x = data$x, y = stem_top + grid::unit(0.7, "mm"), r = grid::unit(0.7, "mm"), gp = loose_gp)
    ), name = "loose_ends")
  }
)

#' @name loose_end_copy_number
#' @title Copy number at each loose end
#'
#' @description
#' Adds a `cn` column with the copy number of the segment containing each loose end,
#' for placing [geom_loose_end()] on a copy number track.
#'
#' @param loose_ends Data frame of loose ends with `seqnames` and `start` columns.
#' @inheritParams junction_copy_number
#'
#' @returns `loose_ends` as a data frame with a `cn` column added
#' @export
#' @keywords data
loose_end_copy_number <- function(loose_ends, segments, cn_column = "total_cn", cn_cap = NULL) {
  loose_ends <- as.data.frame(as_genomic_df(loose_ends))
  loose_ends$cn <- lookup_copy_number(loose_ends$seqnames, loose_ends$start, segments,
                                      cn_column = cn_column, cn_cap = cn_cap)
  loose_ends
}
