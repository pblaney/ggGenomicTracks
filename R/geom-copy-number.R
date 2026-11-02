#' @name geom_copy_number
#' @title Copy number segment track
#'
#' @description
#' Draws copy number segments as thin horizontal bars at their integer copy number, as
#' in gTrack. `xmin`, `xmax`, and `y` default to the `start`, `end`, and `total_cn`
#' columns, as in JaBbA copy number profiles; map `y` to another column such as
#' `major_allele` to change it. Segments are filled by their integer copy number with
#' [scale_fill_copynumber()] unless `fill` is mapped or set, e.g. by state with
#' `aes(fill = cn_state)` and [scale_fill_cn_state()]. Bars use `fill` so `colour` stays
#' free for [geom_junction()] in the same plot.
#'
#' Segments above `cn_cap` are drawn at the cap and marked with a triangle, so a few
#' highly amplified segments do not flatten the rest of the profile. Adjacent segments
#' are joined by thin gray reference connectors, as in JaBbA genome graphs.
#'
#' With `alleles`, the major and minor allele copy numbers are drawn as two
#' semi-transparent series, red and blue, so where they are equal the bars overlap in
#' purple.
#'
#' Also adds a y-axis scale with whole-number breaks starting at zero and a
#' "Copy Number" y-axis title. Add your own `scale_y_continuous()` or `labs()` after it
#' to override them.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()], e.g. `fill`.
#' @param data GRanges object or data frame of copy number segments with `seqnames`,
#'   `start`, `end`, and `total_cn` columns, e.g. a JaBbA `cnv.bed` profile.
#' @param cn_cap Highest copy number shown, default: no cap. Segments above it are drawn
#'   at the cap with a triangle marker.
#' @param bar_height Height of the bars in mm, default: 1.5
#' @param connectors Join adjacent segments with thin gray lines, default: TRUE.
#'   Segments are assumed to come from one sample per layer.
#' @param alleles Draw allele-specific copy number instead of total copy number, by
#'   giving the major and minor allele columns, e.g. `c("major_allele",
#'   "minor_allele")`, or `TRUE` for those JaBbA column names. Default: `NULL` for total
#'   copy number.
#' @param show_first_y_label Label the lowest y-axis break, default: TRUE. See
#'   [label_drop_first()].
#' @param ... Other arguments passed to [ggplot2::layer()], e.g. `fill`, `colour` (bar
#'   outline), or `alpha`.
#' @param na.rm,show.legend,inherit.aes See [ggplot2::layer()].
#'
#' @returns List of a ggplot2 layer, a y-axis scale, a y-axis title, and a copy number
#'   fill scale unless `fill` is given, added to a plot with `+`
#' @export
#' @keywords geom
#'
#' @examples
#' \dontrun{
#' cn <- data.table::fread("sample.jabba.phased.cnv.bed")
#' data.table::setnames(cn, "chrom", "seqnames")
#'
#' ggplot2::ggplot(cn) +
#'   geom_copy_number(cn_cap = 8) +
#'   coord_genome() +
#'   theme_genomic()
#' }
geom_copy_number <- function(mapping = NULL, data = NULL, cn_cap = NULL, bar_height = 1.5,
                             connectors = TRUE, alleles = NULL, show_first_y_label = TRUE, ...,
                             na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  # Default segment aesthetics, overridden by user mapping
  cn_mapping <- ggplot2::aes(xmin = start, xmax = end, y = total_cn)
  for(aes_name in names(mapping)) {
    cn_mapping[[aes_name]] <- mapping[[aes_name]]
  }

  # Allele-specific copy number as major and minor series, each with its own connectors
  layer_params <- list(...)
  if(isTRUE(alleles)) {
    alleles <- c("major_allele", "minor_allele")
  }
  fill_scale <- NULL
  if(!is.null(alleles)) {
    data <- if(is.null(data)) function(plot_data) allele_long(plot_data, alleles) else allele_long(data, alleles)
    if(!"y" %in% names(mapping)) {
      cn_mapping$y <- ggplot2::aes(y = allele_cn)$y
    }
    if(!"fill" %in% c(names(mapping), names(layer_params))) {
      cn_mapping$fill <- ggplot2::aes(fill = allele)$fill
      fill_scale <- ggplot2::scale_fill_manual(name = "Allele", values = c(Major = "#E41A1C", Minor = "#377EB8"))
    }
    if(!"alpha" %in% names(layer_params)) {
      layer_params$alpha <- 0.6
    }
  }

  # Fill by integer copy number, after capping, unless fill is given
  if(is.null(alleles) && !"fill" %in% c(names(cn_mapping), names(layer_params))) {
    cn_mapping$fill <- ggplot2::aes(fill = ggplot2::after_stat(factor(y)))$fill
    fill_scale <- scale_fill_copynumber()
  }

  cn_labels <- scales::label_number(accuracy = 1)
  if(!show_first_y_label) {
    cn_labels <- label_drop_first(cn_labels)
  }

  cn_layer <- ggplot2::layer(geom = GeomCopyNumber, stat = StatCopyNumber, position = ggplot2::PositionIdentity,
                             data = if(is.function(data)) data else as_genomic_df(data), mapping = cn_mapping,
                             show.legend = show.legend, inherit.aes = inherit.aes,
                             params = c(list(cn_cap = cn_cap, bar_height = bar_height, connectors = connectors,
                                             connect_by_fill = !is.null(alleles), na.rm = na.rm),
                                        layer_params))

  list(cn_layer,
       ggplot2::scale_y_continuous(breaks = cn_breaks, labels = cn_labels, limits = c(0, NA),
                                   expand = ggplot2::expansion(add = c(0.4, 0.6))),
       ggplot2::labs(y = "Copy\nNumber"),
       fill_scale)
}

#' @rdname geom_copy_number
#' @format NULL
#' @usage NULL
#' @export
StatCopyNumber <- ggplot2::ggproto("StatCopyNumber", ggplot2::Stat,
  required_aes = c("xmin", "xmax", "y"),
  extra_params = c("na.rm", "cn_cap"),

  # Segments above the cap are drawn at the cap and flagged for a marker
  compute_layer = function(self, data, params, layout) {
    data$capped <- FALSE
    if(!is.null(params$cn_cap)) {
      data$capped <- !is.na(data$y) & data$y > params$cn_cap
      data$y[data$capped] <- params$cn_cap
    }
    data
  }
)

#' @rdname geom_copy_number
#' @format NULL
#' @usage NULL
#' @export
GeomCopyNumber <- ggplot2::ggproto("GeomCopyNumber", ggplot2::Geom,
  required_aes = c("xmin", "xmax", "y"),
  optional_aes = "capped",
  default_aes = ggplot2::aes(fill = "black", colour = NA, linewidth = 0.2, linetype = 1, alpha = NA),
  draw_key = ggplot2::draw_key_rect,

  draw_panel = function(data, panel_params, coord, bar_height = 1.5, connectors = TRUE,
                        connect_by_fill = FALSE, na.rm = FALSE) {
    bars <- coord$transform(data, panel_params)
    cn_grobs <- grid::gList()

    # Reference connectors between segments that touch, within each allele series and
    # never across the boundary of two windows
    if(connectors && nrow(bars) > 1) {
      series <- if(connect_by_fill) bars$fill else rep("all", nrow(bars))
      window <- if(inherits(coord, "CoordLoci")) findInterval(data$xmin, coord$layout$lin_start) else rep(1, nrow(data))
      segment_order <- order(series, bars$xmin)
      ordered <- bars[segment_order, , drop = FALSE]
      ordered_series <- series[segment_order]
      ordered_window <- window[segment_order]
      n <- nrow(ordered)
      touching <- ordered_series[-1] == ordered_series[-n] & ordered_window[-1] == ordered_window[-n] &
        abs(ordered$xmin[-1] - ordered$xmax[-n]) < 0.002 &
        ordered$y[-1] != ordered$y[-n]
      if(any(touching)) {
        connector_x <- (ordered$xmax[-n][touching] + ordered$xmin[-1][touching]) / 2
        cn_grobs <- grid::gList(cn_grobs,
                                grid::segmentsGrob(x0 = connector_x, x1 = connector_x,
                                                   y0 = ordered$y[-n][touching], y1 = ordered$y[-1][touching],
                                                   gp = grid::gpar(col = "gray70", lwd = 0.5)))
      }
    }

    # Bars keep a fixed height in mm at any copy number scale
    cn_grobs <- grid::gList(cn_grobs, grid::rectGrob(x = grid::unit(bars$xmin, "npc"), y = grid::unit(bars$y, "npc"),
                                           width = grid::unit(bars$xmax - bars$xmin, "npc"),
                                           height = grid::unit(bar_height, "mm"), just = c("left", "centre"),
                                           gp = grid::gpar(fill = ggplot2::alpha(bars$fill, bars$alpha),
                                                           col = bars$colour, lwd = bars$linewidth * ggplot2::.pt,
                                                           lty = bars$linetype)))

    # Triangle at the middle of each capped segment
    capped <- data[!is.na(data$capped) & data$capped, , drop = FALSE]
    if(nrow(capped) > 0) {
      markers <- data.frame(x = (capped$xmin + capped$xmax) / 2, y = capped$y,
                            colour = capped$fill, fill = capped$fill, alpha = capped$alpha,
                            shape = 24, size = 2.5, stroke = 0.5, PANEL = capped$PANEL, group = capped$group)
      cn_grobs <- grid::gList(cn_grobs, ggplot2::GeomPoint$draw_panel(markers, panel_params, coord))
    }
    grid::gTree(children = cn_grobs, name = "copy_number")
  }
)

#' @name scale_colour_copynumber
#' @title Integer copy number color scales
#'
#' @description
#' Colors for integer copy number from [devgru::copynumber_palettier()]: blue for 0
#' and 1, gray for 2, and orange to dark red from 3 up to the highest copy number in
#' the plot. Map a discrete copy number, e.g. `aes(fill = factor(total_cn))` with
#' [geom_bed_track()]; [geom_copy_number()] uses the fill scale by default.
#'
#' @param name Legend title, default: "Copy\\nNumber"
#' @param ... Other arguments passed to [ggplot2::discrete_scale()], e.g. `na.value`.
#'
#' @returns ggplot2 scale
#' @export
#' @keywords scale
scale_colour_copynumber <- function(name = "Copy\nNumber", ...) {
  copynumber_scale(aesthetics = "colour", name = name, ...)
}

#' @rdname scale_colour_copynumber
#' @export
scale_fill_copynumber <- function(name = "Copy\nNumber", ...) {
  copynumber_scale(aesthetics = "fill", name = name, ...)
}

#' @name copynumber_scale
#' @title Build an integer copy number scale
#'
#' @description
#' The palette depends on the copy numbers present, so it is built from the scale
#' limits when values are mapped rather than from the number of levels. For internal
#' use only.
#'
#' @param aesthetics Aesthetic of the scale
#' @inheritParams scale_colour_copynumber
#' @keywords internal
copynumber_scale <- function(aesthetics, name, ..., na.value = "black") {
  base_scale <- ggplot2::discrete_scale(aesthetics = aesthetics, palette = function(n) rep(NA_character_, n),
                                        name = name, na.value = na.value, ...)
  ggplot2::ggproto(NULL, base_scale,
    map = function(self, x, limits = self$get_limits()) {
      limits <- limits[!is.na(limits)]
      cn_palette <- devgru::copynumber_palettier(as.numeric(as.character(limits)))
      cn_colors <- stats::setNames(cn_palette$color, as.character(cn_palette$copynumber))
      mapped <- unname(cn_colors[as.character(x)])
      mapped[is.na(mapped)] <- self$na.value
      mapped
    }
  )
}

#' @name allele_long
#' @title Allele-specific copy number in long format
#'
#' @description
#' Stacks the major and minor allele columns into one `allele_cn` column with an
#' `allele` column of `"Major"` or `"Minor"`. For internal use only.
#'
#' @param data Data frame or GRanges of copy number segments
#' @param alleles Names of the major and minor allele columns
#' @keywords internal
allele_long <- function(data, alleles) {
  data <- as.data.frame(as_genomic_df(data))
  missing_columns <- setdiff(alleles, names(data))
  if(length(missing_columns) > 0) {
    cli::cli_abort("Allele column{?s} {.field {missing_columns}} not found in the copy number data.")
  }
  major <- data
  major$allele <- "Major"
  major$allele_cn <- data[[alleles[1]]]
  minor <- data
  minor$allele <- "Minor"
  minor$allele_cn <- data[[alleles[2]]]
  long_data <- rbind(major, minor)
  long_data$allele <- factor(long_data$allele, levels = c("Major", "Minor"))
  long_data
}

#' @name cn_breaks
#' @title Copy number axis breaks
#'
#' @description
#' Every whole number for small copy number ranges, otherwise up to five pretty whole
#' numbers. For internal use only.
#'
#' @param limits Axis limits
#' @keywords internal
cn_breaks <- function(limits) {
  if(diff(limits) <= 6) {
    return(seq(ceiling(limits[1]), floor(limits[2])))
  }
  breaks <- scales::breaks_pretty(n = 5)(limits)
  breaks[breaks == round(breaks)]
}

#' @name scale_colour_cn_state
#' @title Copy number state color scales
#'
#' @description
#' Colors for the JaBbA copy number states, from homozygous deletion in dark blue to
#' amplification in dark red, with neutral in gray. Use the fill scale with
#' [geom_copy_number()] or with [geom_bed_track()] to show states across many samples.
#'
#' @param name Legend title, default: "CN State"
#' @param ... Other arguments passed to [ggplot2::scale_colour_manual()] or
#'   [ggplot2::scale_fill_manual()], e.g. `values` to change colors.
#'
#' @returns ggplot2 scale
#' @export
#' @keywords scale
scale_colour_cn_state <- function(name = "CN State", ...) {
  do.call(ggplot2::scale_colour_manual, c(list(name = name, ...), cn_state_scale_defaults(...)))
}

#' @rdname scale_colour_cn_state
#' @export
scale_fill_cn_state <- function(name = "CN State", ...) {
  do.call(ggplot2::scale_fill_manual, c(list(name = name, ...), cn_state_scale_defaults(...)))
}

#' @name cn_state_scale_defaults
#' @title Default arguments of the copy number state scales
#'
#' @description
#' State colors and legend order, unless given in `...`. For internal use only.
#'
#' @param ... Arguments given to the scale
#' @keywords internal
cn_state_scale_defaults <- function(...) {
  given <- names(list(...))
  defaults <- list(values = cn_state_colors, breaks = names(cn_state_colors))
  defaults[setdiff(names(defaults), given)]
}

# JaBbA copy number states, from deletion to amplification
cn_state_colors <- c("homozygous_deletion" = "#08306B", "loss" = "#4A90C8", "neutral" = "gray60",
                     "gain" = "#F4845F", "amplification" = "#A50F15")
