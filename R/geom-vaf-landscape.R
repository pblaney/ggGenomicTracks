#' @name geom_vaf_landscape
#' @title Variant allele frequency track
#'
#' @description
#' Draws each mutation as a point at its position and variant allele frequency (VAF),
#' colored by the local density of points so clusters of mutations at similar VAF
#' stand out. `x` and `y` default to the `start` and `VAF` columns. Density is
#' estimated from the mutations within the [coord_locus()] window.
#'
#' Also adds a "Density" color scale labeled from low to high, a y-axis fixed to 0-1,
#' and a "VAF" y-axis title. Add your own scales or `labs()` after it to override them.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()].
#' @param data GRanges object or data frame of mutations with `seqnames`, `start`, and
#'   `VAF` columns, e.g. a MAF with `VAF = t_alt_count / t_depth`.
#' @param point_size Size of the points in mm, default: 3
#' @param background Fill color of the panel behind the points, default: `"gray50"`,
#'   which gives the light end of the density colors contrast. Use `NA` for none.
#' @param method,adjust Density estimation method and bandwidth adjustment, see
#'   [ggpointdensity::stat_pointdensity()]. Default: `"kde2d"`, 1
#' @param show_first_y_label Label the lowest y-axis break, default: TRUE. See
#'   [label_drop_first()].
#' @param ... Other arguments passed to [ggpointdensity::stat_pointdensity()], e.g.
#'   `alpha` or `shape`.
#' @param na.rm,show.legend,inherit.aes See [ggplot2::layer()].
#'
#' @returns List of ggplot2 layers, scales, and a y-axis title, added to a plot with `+`
#' @export
#' @keywords geom
#'
#' @examples
#' \dontrun{
#' snvs$VAF <- snvs$t_alt_count / snvs$t_depth
#'
#' ggplot2::ggplot(snvs) +
#'   geom_vaf_landscape() +
#'   coord_locus("chr1:203296001-203313000") +
#'   theme_genomic()
#' }
geom_vaf_landscape <- function(mapping = NULL, data = NULL, point_size = 3, background = "gray50",
                               method = "kde2d", adjust = 1, show_first_y_label = TRUE, ...,
                               na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  # Default position aesthetics, overridden by user mapping
  vaf_mapping <- ggplot2::aes(x = start, y = VAF)
  for(aes_name in names(mapping)) {
    vaf_mapping[[aes_name]] <- mapping[[aes_name]]
  }

  vaf_labels <- scales::label_number(accuracy = 0.1)
  if(!show_first_y_label) {
    vaf_labels <- label_drop_first(vaf_labels)
  }

  # Background drawn as a layer rather than a theme setting so theme_genomic() keeps it
  background_layer <- NULL
  if(!is.na(background)) {
    background_layer <- ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                                          fill = background)
  }

  list(background_layer,
       ggpointdensity::stat_pointdensity(mapping = vaf_mapping, data = as_genomic_df(data),
                                         method = method, adjust = adjust, size = point_size, ...,
                                         na.rm = na.rm, show.legend = show.legend,
                                         inherit.aes = inherit.aes),
       # Raw density values depend on the window size, so only the ends are labeled
       ggplot2::scale_colour_gradientn(name = "Density",
                                       colours = paletteer::paletteer_c("grDevices::Lajolla", n = 15),
                                       breaks = function(limits) limits, labels = c("Low", "High")),
       ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), labels = vaf_labels,
                                   expand = c(0.01, 0.01)),
       ggplot2::labs(y = "VAF"))
}
