#' @name geom_mutation_count
#' @title Mutation count lollipop track
#'
#' @description
#' Counts mutations at each position and draws the counts as lollipops: a stem from
#' zero and a point sized by the count. `x` defaults to the `start` column. Counts are
#' made within each group, so mapping `fill` or `group` (e.g. to `Reference_Allele` or
#' a sample cohort) counts each level separately.
#'
#' Also adds a y-axis scale with whole-number breaks starting at zero and a
#' "Mutation Count" y-axis title. Add your own `scale_y_continuous()` or `labs()` after
#' it to override them.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()], e.g. `fill`.
#' @param data GRanges object or data frame of mutations with `seqnames` and `start`
#'   columns, one row per mutation, e.g. from a MAF or VCF read with devgru.
#' @param size_range Range of point sizes in mm for the smallest and largest counts,
#'   default: `c(1.5, 4)`. Point area grows with the count. Map `size` to
#'   `after_stat(count)` instead to get a size legend.
#' @param stem_colour,stem_alpha,stem_linewidth Color, opacity, and width in mm of the
#'   stems.
#' @param ... Other arguments passed to [ggplot2::layer()], e.g. `fill`, `colour`
#'   (point outline), `shape`, or `alpha`.
#' @param na.rm,show.legend,inherit.aes See [ggplot2::layer()].
#'
#' @section Computed variables:
#' * `count`: number of mutations at the position within the group.
#'
#' @returns List of a ggplot2 layer, a y-axis scale, and a y-axis title, added to a
#'   plot with `+`
#' @export
#' @keywords geom
#'
#' @examples
#' \dontrun{
#' ggplot2::ggplot(snvs) +
#'   geom_mutation_count(fill = "dodgerblue") +
#'   coord_locus("chr3:187734001-187751500", flank = 5000) +
#'   theme_genomic()
#'
#' # Count each cohort separately
#' ggplot2::ggplot(snvs) +
#'   geom_mutation_count(ggplot2::aes(fill = cohort), alpha = 0.8) +
#'   coord_locus("chr3:187734001-187751500", flank = 5000) +
#'   theme_genomic()
#' }
geom_mutation_count <- function(mapping = NULL, data = NULL, size_range = c(1.5, 4),
                                stem_colour = "gray10", stem_alpha = 0.65, stem_linewidth = 0.5,
                                ..., na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  # Default position aesthetic, overridden by user mapping
  count_mapping <- ggplot2::aes(x = start)
  for(aes_name in names(mapping)) {
    count_mapping[[aes_name]] <- mapping[[aes_name]]
  }

  count_layer <- ggplot2::layer(geom = GeomMutationCount, stat = StatMutationCount,
                                position = ggplot2::PositionIdentity,
                                data = as_genomic_df(data), mapping = count_mapping,
                                show.legend = show.legend, inherit.aes = inherit.aes,
                                params = list(size_range = size_range, stem_colour = stem_colour,
                                              stem_alpha = stem_alpha, stem_linewidth = stem_linewidth,
                                              na.rm = na.rm, ...))

  list(count_layer,
       ggplot2::scale_y_continuous(breaks = integer_breaks, labels = scales::label_number(accuracy = 1),
                                   expand = ggplot2::expansion(mult = c(0, 0.15))),
       ggplot2::labs(y = "Mutation\nCount"))
}

#' @rdname geom_mutation_count
#' @format NULL
#' @usage NULL
#' @export
StatMutationCount <- ggplot2::ggproto("StatMutationCount", ggplot2::Stat,
  required_aes = "x",
  default_aes = ggplot2::aes(y = ggplot2::after_stat(count)),

  setup_params = function(data, params) {
    # Point sizes are scaled by the largest count across the whole layer
    params$max_count <- max(table(paste(data$PANEL, data$group, data$x)))
    params
  },

  compute_group = function(data, scales, size_range = c(1.5, 4), max_count = 1) {
    position_index <- match(data$x, unique(data$x))
    counts <- data[!duplicated(position_index), , drop = FALSE]
    counts$count <- tabulate(position_index)
    counts$ymin <- 0

    # Point area grows with count, from the smallest size at one mutation
    if(!"size" %in% names(counts)) {
      counts$size <- scales::rescale(sqrt(counts$count), to = size_range,
                                     from = c(1, sqrt(max(max_count, 2))))
    }
    counts
  }
)

#' @rdname geom_mutation_count
#' @format NULL
#' @usage NULL
#' @export
GeomMutationCount <- ggplot2::ggproto("GeomMutationCount", ggplot2::Geom,
  required_aes = c("x", "y"),
  optional_aes = "ymin",
  default_aes = ggplot2::aes(shape = 21, colour = NA, fill = "gray70", size = 2,
                             alpha = NA, stroke = 0.5),
  draw_key = ggplot2::draw_key_point,

  draw_panel = function(data, panel_params, coord, stem_colour = "gray10", stem_alpha = 0.65,
                        stem_linewidth = 0.5, na.rm = FALSE) {
    stems <- data.frame(x = data$x, xend = data$x, y = data$ymin, yend = data$y,
                        colour = stem_colour, alpha = stem_alpha, linewidth = stem_linewidth,
                        linetype = 1, PANEL = data$PANEL, group = data$group)

    # Smaller points are drawn last so they stay visible where positions overlap
    points <- data[order(-data$size), , drop = FALSE]

    grid::gTree(children = grid::gList(
      ggplot2::GeomSegment$draw_panel(stems, panel_params, coord, lineend = "butt"),
      ggplot2::GeomPoint$draw_panel(points, panel_params, coord, na.rm = na.rm)
    ), name = "mutation_count")
  }
)

#' @name integer_breaks
#' @title Whole-number axis breaks
#'
#' @description
#' Up to three pretty breaks, keeping whole numbers only. For internal use only.
#'
#' @param limits Axis limits
#' @keywords internal
integer_breaks <- function(limits) {
  breaks <- scales::breaks_pretty(n = 3)(limits)
  breaks[breaks == round(breaks)]
}
