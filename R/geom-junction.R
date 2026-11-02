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
#' @param arc_shape Shape of the arcs. `"arch"` (default) draws symmetric arches;
#'   `"strand"` draws arcs leaving each breakpoint toward the side it joins, from the
#'   `strand1` and `strand2` columns, as in gTrack: deletion-like junctions arch inward,
#'   duplication-like junctions loop outward, and inversions hook to one side, so short
#'   fold-back inversions show as hooks. Stubs lean the same way.
#' @param breakpoint_marks Mark each breakpoint with a short tick pointing to the side
#'   it joins, from the `strand1` and `strand2` columns, as in PCAWG figures. Default:
#'   FALSE
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
                          scale_by_span = TRUE, arc_shape = c("arch", "strand"), partner_label_size = 2.5,
                          breakpoint_marks = FALSE, ..., na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  arc_direction <- match.arg(arc_direction)
  arc_shape <- match.arg(arc_shape)

  # Default breakpoint aesthetics, with the chromosome of an out of view partner as the
  # stub label, overridden by user mapping
  junction_mapping <- ggplot2::aes(x = start1, xend = start2,
                                   label = ifelse(is.na(start2), as.character(seqnames2), as.character(seqnames1)))
  if(breakpoint_marks || arc_shape == "strand") {
    junction_mapping$strand <- ggplot2::aes(strand = strand1)$strand
    junction_mapping$strandend <- ggplot2::aes(strandend = strand2)$strandend
  }
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
                                                 scale_by_span = scale_by_span, arc_shape = arc_shape,
                                                 partner_label_size = partner_label_size,
                                                 breakpoint_marks = breakpoint_marks, na.rm = na.rm, ...))
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
  optional_aes = c("label", "strand", "strandend"),
  default_aes = ggplot2::aes(colour = "black", linewidth = 0.5, linetype = 1, alpha = NA),
  draw_key = ggplot2::draw_key_path,

  # A missing breakpoint marks a junction partly out of view, drawn as a stub
  handle_na = function(self, data, params) {
    data[!(is.na(data$x) & is.na(data$xend)), , drop = FALSE]
  },

  draw_panel = function(data, panel_params, coord, arc_height = 0.9, arc_direction = "up",
                        scale_by_span = TRUE, arc_shape = "arch", partner_label_size = 2.5,
                        breakpoint_marks = FALSE, na.rm = FALSE) {
    data <- coord$transform(data, panel_params)
    direction <- if(arc_direction == "up") 1 else -1
    line_gp <- function(rows) {
      grid::gpar(col = ggplot2::alpha(rows$colour, rows$alpha), lwd = rows$linewidth * ggplot2::.pt,
                 lty = rows$linetype, lineend = "butt")
    }
    junction_grobs <- grid::gList()

    # Direction out of the joined segment: right from a + breakpoint, which joins the
    # segment ending at it, left from a - breakpoint, or 0 for arches
    strand_side <- function(strand) {
      if(arc_shape != "strand" || is.null(strand)) {
        return(rep(0, nrow(data)))
      }
      side <- c("+" = 1, "-" = -1)[as.character(strand)]
      ifelse(is.na(side), 0, side)
    }
    data$side <- strand_side(data$strand)
    data$sideend <- strand_side(data$strandend)

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

      # Arches have control points above the breakpoints; strand arcs shift them toward
      # the joined side, continuing out of the segment end
      reach <- pmax(0.3 * span, 0.01)
      control_x <- arcs$x + arcs$side * reach
      control_xend <- arcs$xend + arcs$sideend * reach

      t <- seq(0, 1, length.out = 40)
      b0 <- (1 - t)^3
      b1 <- 3 * (1 - t)^2 * t
      b2 <- 3 * (1 - t) * t^2
      b3 <- t^3
      arc_x <- outer(arcs$x, b0) + outer(control_x, b1) + outer(control_xend, b2) + outer(arcs$xend, b3)
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
      stub_side <- ifelse(is.na(stubs$x), stubs$sideend, stubs$side)
      stub_top <- stub_y + direction * arc_height * 0.5
      stub_tip <- grid::unit(stub_x, "npc") + grid::unit(stub_side * 2, "mm")
      junction_grobs <- grid::gList(junction_grobs,
                                    grid::segmentsGrob(x0 = stub_x, x1 = stub_tip, y0 = stub_y, y1 = stub_top,
                                                       arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "open"),
                                                       gp = line_gp(stubs)))
      if(partner_label_size > 0 && "label" %in% names(stubs)) {
        junction_grobs <- grid::gList(junction_grobs,
                                      grid::textGrob(label = stubs$label, x = stub_tip,
                                                     y = grid::unit(stub_top, "npc") + grid::unit(direction * 0.5, "mm"),
                                                     vjust = if(direction == 1) 0 else 1, rot = 0,
                                                     gp = grid::gpar(col = ggplot2::alpha(stubs$colour, stubs$alpha),
                                                                     fontsize = partner_label_size * ggplot2::.pt)))
      }
    }
    # Ticks at each breakpoint in view, pointing left for + and right for - breakpoints
    if(breakpoint_marks && all(c("strand", "strandend") %in% names(data))) {
      ends <- rbind(data.frame(x = data$x, y = data$y, strand = data$strand, colour = data$colour,
                               alpha = data$alpha, linewidth = data$linewidth),
                    data.frame(x = data$xend, y = data$yend, strand = data$strandend, colour = data$colour,
                               alpha = data$alpha, linewidth = data$linewidth))
      ends <- ends[!is.na(ends$x) & ends$strand %in% c("+", "-"), , drop = FALSE]
      if(nrow(ends) > 0) {
        tick_side <- ifelse(ends$strand == "+", -1, 1)
        junction_grobs <- grid::gList(junction_grobs,
                                      grid::segmentsGrob(x0 = grid::unit(ends$x, "npc"),
                                                         x1 = grid::unit(ends$x, "npc") + grid::unit(tick_side * 1.5, "mm"),
                                                         y0 = ends$y, y1 = ends$y,
                                                         gp = grid::gpar(col = ggplot2::alpha(ends$colour, ends$alpha),
                                                                         lwd = ends$linewidth * 3 * ggplot2::.pt,
                                                                         lineend = "butt")))
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
  for(anchor in c("1", "2")) {
    position <- junctions[[paste0("start", anchor)]]
    strand_column <- paste0("strand", anchor)
    if(strand_column %in% names(junctions)) {
      position <- ifelse(junctions[[strand_column]] == "-", position + 1, position)
    }
    junctions[[paste0("cn", anchor)]] <- lookup_copy_number(junctions[[paste0("seqnames", anchor)]], position,
                                                            segments, cn_column = cn_column, cn_cap = cn_cap)
  }
  junctions
}

#' @name lookup_copy_number
#' @title Copy number of the segment containing each position
#'
#' @description
#' For internal use only.
#'
#' @param seqnames,position Chromosomes and positions to look up
#' @inheritParams junction_copy_number
#' @returns Copy number per position, `NA` where no segment contains it
#' @keywords internal
lookup_copy_number <- function(seqnames, position, segments, cn_column = "total_cn", cn_cap = NULL) {
  segments <- as.data.frame(as_genomic_df(segments))
  segment_key <- chrom_keys(segments$seqnames)
  position_key <- chrom_keys(seqnames)

  position_cn <- rep(NA_real_, length(position))
  for(chrom in intersect(unique(position_key), segment_key)) {
    chrom_segments <- segments[segment_key == chrom, , drop = FALSE]
    chrom_segments <- chrom_segments[order(chrom_segments$start), , drop = FALSE]
    on_chrom <- which(position_key == chrom)
    segment_index <- findInterval(position[on_chrom], chrom_segments$start)
    found <- segment_index > 0
    found[found] <- position[on_chrom][found] <= chrom_segments$end[segment_index[found]]
    position_cn[on_chrom[found]] <- chrom_segments[[cn_column]][segment_index[found]]
  }
  if(!is.null(cn_cap)) {
    position_cn <- pmin(position_cn, cn_cap)
  }
  position_cn
}

#' @name annotate_junctions
#' @title Classify junctions by orientation and complex event
#'
#' @description
#' Adds two columns for coloring [geom_junction()]:
#' * `sv_class`: orientation from the breakpoint strands, `"DEL-like"` (+/-),
#'   `"DUP-like"` (-/+), `"INV-h2h"` (+/+), or `"INV-t2t"` (-/-), for both cis and trans
#'   junctions. Use with [scale_colour_sv_class()].
#' * `event`: the JaBbA complex event the junction belongs to, from the event columns
#'   (`chromothripsis`, `chromoplexy`, `templated_insertion`, `rigma`, `pyrgo`, `tyfonas`,
#'   `ecdna`, `bfb`), or `"Simple"`; and `event_id` with the event label, e.g. `"CPX-4"`.
#'   Use with [scale_colour_jabba_event()].
#'
#' @param junctions Data frame of junctions with `strand1` and `strand2` columns and
#'   optionally JaBbA event columns, e.g. from `devgru::read_bedpe_file(output_type =
#'   "dt")` of a JaBbA `junctions.bedpe`.
#'
#' @returns `junctions` as a data frame with `sv_class`, `event`, and `event_id` added
#' @export
#' @keywords data
annotate_junctions <- function(junctions) {
  junctions <- as.data.frame(junctions)
  if(all(c("strand1", "strand2") %in% names(junctions))) {
    orientation <- paste0(junctions$strand1, junctions$strand2)
    junctions$sv_class <- factor(c("+-" = "DEL-like", "-+" = "DUP-like", "++" = "INV-h2h", "--" = "INV-t2t")[orientation],
                                 levels = names(sv_class_colors$jabba))
  }

  # First event column with a label, in order of the JaBbA event columns
  event_columns <- intersect(names(jabba_event_labels), names(junctions))
  junctions$event <- "Simple"
  junctions$event_id <- NA_character_
  for(column in rev(event_columns)) {
    labeled <- !is.na(junctions[[column]]) & as.character(junctions[[column]]) != ""
    junctions$event[labeled] <- jabba_event_labels[[column]]
    junctions$event_id[labeled] <- as.character(junctions[[column]][labeled])
  }
  junctions$event <- factor(junctions$event, levels = names(jabba_event_colors))
  junctions
}

#' @name scale_colour_sv_class
#' @title Junction orientation color scale
#'
#' @description
#' Colors for junction orientation classes from [annotate_junctions()]. The `"jabba"`
#' palette follows the JaBbA paper, with DEL-like junctions in cyan and DUP-like in
#' magenta; the `"pcawg"` palette follows the PCAWG structural variant paper, with
#' deletion-type in blue, duplication-type in red, and head-to-head and tail-to-tail
#' inversions in green and purple.
#'
#' @param palette `"jabba"` (default) or `"pcawg"`
#' @param name Legend title, default: "Junction"
#' @param ... Other arguments passed to [ggplot2::scale_colour_manual()].
#'
#' @returns ggplot2 scale
#' @export
#' @keywords scale
scale_colour_sv_class <- function(palette = c("jabba", "pcawg"), name = "Junction", ...) {
  palette <- match.arg(palette)
  ggplot2::scale_colour_manual(name = name, values = sv_class_colors[[palette]],
                               breaks = names(sv_class_colors[[palette]]), ...)
}

#' @name scale_colour_jabba_event
#' @title JaBbA complex event color scale
#'
#' @description
#' Colors for the complex event classes from [annotate_junctions()], with junctions
#' outside complex events in light gray so events stand out.
#'
#' @param name Legend title, default: "Event"
#' @param ... Other arguments passed to [ggplot2::scale_colour_manual()].
#'
#' @returns ggplot2 scale
#' @export
#' @keywords scale
scale_colour_jabba_event <- function(name = "Event", ...) {
  ggplot2::scale_colour_manual(name = name, values = jabba_event_colors, ...)
}

# JaBbA event columns and their labels
jabba_event_labels <- c(chromothripsis = "Chromothripsis", chromoplexy = "Chromoplexy",
                        templated_insertion = "TIC", rigma = "Rigma", pyrgo = "Pyrgo",
                        tyfonas = "Tyfonas", ecdna = "ecDNA", bfb = "BFB")

jabba_event_colors <- c("Chromothripsis" = "#7B3294", "Chromoplexy" = "#E66101", "TIC" = "#1B9E77",
                        "Rigma" = "#D01C8B", "Pyrgo" = "#B8860B", "Tyfonas" = "#E41A1C",
                        "ecDNA" = "#A6611A", "BFB" = "#2166AC", "Simple" = "gray75")

sv_class_colors <- list(
  jabba = c("DEL-like" = "#00B0F0", "DUP-like" = "#E000E0", "INV-h2h" = "#2CA25F", "INV-t2t" = "#7B3294"),
  pcawg = c("DEL-like" = "#2166AC", "DUP-like" = "#D7301F", "INV-h2h" = "#1A9850", "INV-t2t" = "#762A83")
)

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
