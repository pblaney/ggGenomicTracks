#' @name geom_junction
#' @title Structural variant junction track
#'
#' @description
#' Draws structural variant junctions as arcs between their two breakpoints. Pair with
#' [coord_loci()] to connect breakpoints in different windows, including trans
#' junctions between chromosomes, or [coord_genome()] for a genome-wide view. A
#' junction with one breakpoint out of view is drawn as a stub at the breakpoint in
#' view, labeled with the chromosome of its partner.
#'
#' `x` and `xend` default to the `start1` and `start2` columns of BEDPE data, as from
#' `devgru::read_bedpe_file(output_type = "dt")`. On its own, the track draws arcs up
#' from the bottom of the panel, or down from the top, taller for junctions spanning
#' more of the axis. To
#' draw junctions on a copy number track, map `y` and `yend` to the copy number at each
#' breakpoint, e.g. from [junction_copy_number()], and arcs connect those heights.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()], e.g.
#'   `colour = type`, or `y` and `yend` to anchor arcs at copy number.
#' @param data Data frame of junctions with `seqnames1`, `start1`, `seqnames2`, and
#'   `start2` columns.
#' @param arc_height Height of the tallest arcs as a fraction of the panel height.
#'   Default: 0.9 on their own, 0.25 when `y` is mapped.
#' @param arc_direction Draw arcs `"up"` (default) or `"down"` from their breakpoints.
#' @param scale_by_span Make arcs spanning more of the axis taller, default: TRUE
#' @param partner_label_size Size of the partner chromosome labels on stubs in mm, or 0
#'   for none. Default: 2.5
#' @param ... Other arguments passed to [ggplot2::layer()], e.g. `colour`, `linewidth`,
#'   or `alpha`.
#' @param na.rm,show.legend,inherit.aes See [ggplot2::layer()].
#'
#' @returns List of a ggplot2 layer and, unless `y` is mapped, a y-axis scale and title,
#'   added to a plot with `+`
#' @export
#' @keywords geom
#'
#' @examples
#' \dontrun{
#' junctions <- devgru::read_bedpe_file("sample.jabba.junctions.bedpe", output_type = "dt")
#' windows <- c(MYC = "chr8:126000001-130000000", IGH = "chr14:105000001-107000000")
#'
#' ggplot2::ggplot(junctions) +
#'   geom_junction(ggplot2::aes(colour = type)) +
#'   scale_colour_sv_type() +
#'   coord_loci(windows) +
#'   theme_genomic()
#'
#' # Junctions joining copy number segments
#' ggplot2::ggplot() +
#'   geom_copy_number(data = cn, cn_cap = 8) +
#'   geom_junction(data = junction_copy_number(junctions, cn, cn_cap = 8),
#'                 ggplot2::aes(y = cn1, yend = cn2, colour = type)) +
#'   scale_colour_sv_type() +
#'   coord_loci(windows) +
#'   theme_genomic()
#' }
geom_junction <- function(mapping = NULL, data = NULL, arc_height = NULL, arc_direction = c("up", "down"),
                          scale_by_span = TRUE, partner_label_size = 2.5, ...,
                          na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  arc_direction <- match.arg(arc_direction)

  # Default breakpoint aesthetics, with the chromosome of an out of view partner as the
  # stub label, overridden by user mapping
  junction_mapping <- ggplot2::aes(x = start1, xend = start2,
                                   label = ifelse(is.na(start2), as.character(seqnames2), as.character(seqnames1)))
  for(aes_name in names(mapping)) {
    junction_mapping[[aes_name]] <- mapping[[aes_name]]
  }

  # Anchored to copy number when y is mapped, otherwise arcs rise from the panel bottom,
  # or hang from the top when drawn down
  anchored <- "y" %in% names(mapping)
  if(!anchored) {
    junction_mapping$y <- if(arc_direction == "up") ggplot2::aes(y = 0)$y else ggplot2::aes(y = 1)$y
  }
  if(!"yend" %in% names(junction_mapping)) {
    junction_mapping$yend <- junction_mapping$y
  }
  if(is.null(arc_height)) {
    arc_height <- if(anchored) 0.25 else 0.9
  }

  junction_layer <- ggplot2::layer(geom = GeomJunction, stat = ggplot2::StatIdentity,
                                   position = ggplot2::PositionIdentity,
                                   data = data, mapping = junction_mapping,
                                   show.legend = show.legend, inherit.aes = inherit.aes,
                                   params = list(arc_height = arc_height, arc_direction = arc_direction,
                                                 scale_by_span = scale_by_span,
                                                 partner_label_size = partner_label_size, na.rm = na.rm, ...))
  if(anchored) {
    return(junction_layer)
  }

  list(junction_layer,
       ggplot2::scale_y_continuous(limits = c(0, 1), breaks = NULL, expand = c(0, 0)),
       ggplot2::labs(y = "SVs"))
}

#' @rdname geom_junction
#' @format NULL
#' @usage NULL
#' @export
GeomJunction <- ggplot2::ggproto("GeomJunction", ggplot2::Geom,
  required_aes = c("x", "xend", "y", "yend"),
  optional_aes = "label",
  default_aes = ggplot2::aes(colour = "black", linewidth = 0.5, linetype = 1, alpha = NA),
  draw_key = ggplot2::draw_key_path,

  # A missing breakpoint marks a junction partly out of view, drawn as a stub
  handle_na = function(self, data, params) {
    data[!(is.na(data$x) & is.na(data$xend)), , drop = FALSE]
  },

  draw_panel = function(data, panel_params, coord, arc_height = 0.9, arc_direction = "up",
                        scale_by_span = TRUE, partner_label_size = 2.5, na.rm = FALSE) {
    data <- coord$transform(data, panel_params)
    direction <- if(arc_direction == "up") 1 else -1
    line_gp <- function(rows) {
      grid::gpar(col = ggplot2::alpha(rows$colour, rows$alpha), lwd = rows$linewidth * ggplot2::.pt,
                 lty = rows$linetype, lineend = "butt")
    }
    junction_grobs <- grid::gList()

    # Arcs between breakpoints in view, as cubic curves peaking at the arc height
    arcs <- data[!is.na(data$x) & !is.na(data$xend), , drop = FALSE]
    if(nrow(arcs) > 0) {
      span <- pmin(abs(arcs$xend - arcs$x), 1)
      height <- arc_height * if(scale_by_span) 0.3 + 0.7 * sqrt(span) else 1
      base <- if(direction == 1) pmax(arcs$y, arcs$yend) else pmin(arcs$y, arcs$yend)
      # Peaks are kept inside the panel
      peak <- base + direction * height
      peak <- if(direction == 1) pmax(pmin(peak, 0.98), base) else pmin(pmax(peak, 0.02), base)
      control <- (peak - 0.25 * (arcs$y + arcs$yend) / 2) / 0.75

      t <- seq(0, 1, length.out = 40)
      b0 <- (1 - t)^3
      b1 <- 3 * (1 - t)^2 * t
      b2 <- 3 * (1 - t) * t^2
      b3 <- t^3
      arc_x <- outer(arcs$x, b0 + b1) + outer(arcs$xend, b2 + b3)
      arc_y <- outer(arcs$y, b0) + outer(control, b1 + b2) + outer(arcs$yend, b3)
      junction_grobs <- grid::gList(junction_grobs,
                                    grid::polylineGrob(x = as.vector(t(arc_x)), y = as.vector(t(arc_y)),
                                                       id = rep(seq_len(nrow(arcs)), each = length(t)),
                                                       gp = line_gp(arcs)))
    }

    # Stubs at breakpoints whose partner is out of view, labeled with its chromosome
    stubs <- data[is.na(data$x) != is.na(data$xend), , drop = FALSE]
    if(nrow(stubs) > 0) {
      stub_x <- ifelse(is.na(stubs$x), stubs$xend, stubs$x)
      stub_y <- ifelse(is.na(stubs$x), stubs$yend, stubs$y)
      stub_top <- stub_y + direction * arc_height * 0.5
      junction_grobs <- grid::gList(junction_grobs,
                                    grid::segmentsGrob(x0 = stub_x, x1 = stub_x, y0 = stub_y, y1 = stub_top,
                                                       arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "open"),
                                                       gp = line_gp(stubs)))
      if(partner_label_size > 0 && "label" %in% names(stubs)) {
        junction_grobs <- grid::gList(junction_grobs,
                                      grid::textGrob(label = stubs$label, x = stub_x,
                                                     y = grid::unit(stub_top, "npc") + grid::unit(direction * 0.5, "mm"),
                                                     vjust = if(direction == 1) 0 else 1, rot = 0,
                                                     gp = grid::gpar(col = ggplot2::alpha(stubs$colour, stubs$alpha),
                                                                     fontsize = partner_label_size * ggplot2::.pt)))
      }
    }
    grid::gTree(children = junction_grobs, name = "junctions")
  }
)

#' @name junction_copy_number
#' @title Copy number at each junction breakpoint
#'
#' @description
#' Adds the copy number on the joined side of each breakpoint, `cn1` and `cn2`, for
#' anchoring [geom_junction()] arcs on a copy number track. A breakpoint on the `+`
#' strand joins the segment ending at it, on its left, and one on the `-` strand joins
#' the segment starting after it, on its right; without strand columns the segment
#' containing the breakpoint is used.
#'
#' @param junctions Data frame of junctions with `seqnames1`, `start1`, `seqnames2`, and
#'   `start2` columns, and optionally `strand1` and `strand2`.
#' @param segments GRanges or data frame of copy number segments with `seqnames`,
#'   `start`, and `end` columns, for a single sample.
#' @param cn_column Column of `segments` with the copy number, default: `"total_cn"`
#' @param cn_cap Highest copy number, matching `cn_cap` of [geom_copy_number()] so arcs
#'   meet capped segments. Default: no cap.
#'
#' @returns `junctions` as a data frame with `cn1` and `cn2` columns added
#' @export
#' @keywords data
junction_copy_number <- function(junctions, segments, cn_column = "total_cn", cn_cap = NULL) {
  junctions <- as.data.frame(junctions)
  segments <- as.data.frame(as_genomic_df(segments))
  segment_key <- stringr::str_remove(as.character(segments$seqnames), "^chr")

  for(anchor in c("1", "2")) {
    position <- junctions[[paste0("start", anchor)]]
    strand_column <- paste0("strand", anchor)
    if(strand_column %in% names(junctions)) {
      position <- ifelse(junctions[[strand_column]] == "-", position + 1, position)
    }
    anchor_key <- stringr::str_remove(as.character(junctions[[paste0("seqnames", anchor)]]), "^chr")

    anchor_cn <- rep(NA_real_, nrow(junctions))
    for(chrom in intersect(unique(anchor_key), segment_key)) {
      chrom_segments <- segments[segment_key == chrom, , drop = FALSE]
      chrom_segments <- chrom_segments[order(chrom_segments$start), , drop = FALSE]
      on_chrom <- which(anchor_key == chrom)
      segment_index <- findInterval(position[on_chrom], chrom_segments$start)
      found <- segment_index > 0
      found[found] <- position[on_chrom][found] <= chrom_segments$end[segment_index[found]]
      anchor_cn[on_chrom[found]] <- chrom_segments[[cn_column]][segment_index[found]]
    }
    if(!is.null(cn_cap)) {
      anchor_cn <- pmin(anchor_cn, cn_cap)
    }
    junctions[[paste0("cn", anchor)]] <- anchor_cn
  }
  junctions
}

#' @name scale_colour_sv_type
#' @title Structural variant type color scale
#'
#' @description
#' Colors for deletions (`DEL`), duplications (`DUP`), inversions (`INV`), and
#' translocations (`TRA`).
#'
#' @param name Legend title, default: "SV Type"
#' @param ... Other arguments passed to [ggplot2::scale_colour_manual()], e.g. `values`
#'   to change colors.
#'
#' @returns ggplot2 scale
#' @export
#' @keywords scale
scale_colour_sv_type <- function(name = "SV Type", ...) {
  scale_args <- list(name = name, ...)
  if(!"values" %in% names(scale_args)) {
    scale_args$values <- sv_type_colors
  }
  if(!"breaks" %in% names(scale_args)) {
    scale_args$breaks <- names(sv_type_colors)
  }
  do.call(ggplot2::scale_colour_manual, scale_args)
}

# JAMA palette colors, as used by the earlier breakpoint track
sv_type_colors <- c("DEL" = "#374E55", "DUP" = "#DF8F44", "INV" = "#00A1D5", "TRA" = "#B24745")
