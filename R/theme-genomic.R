#' @name theme_genomic
#' @title Theme for genomic track plots
#'
#' @description
#' Minimal theme for a single genomic track: blank panel, y-axis line with a
#' horizontal title, and tight margins for stacking tracks. The x-axis style is set
#' with `x_axis`; use `"full"` on the bottom track of a stack and `"ticks"` or
#' `"blank"` on the tracks above it.
#'
#' @param x_axis Style of the x-axis, one of `"full"`, `"ticks"`, or `"blank"`.
#' @param base_size Base font size, default: 11
#'
#' @returns ggplot2 theme object
#' @export
#' @keywords theme
theme_genomic <- function(x_axis = c("full", "ticks", "blank"), base_size = 11) {
  x_axis <- match.arg(x_axis)

  # Base track theme
  track_theme <- ggplot2::theme_grey(base_size = base_size) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   panel.grid = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_line(linewidth = 0.3),
                   axis.title.y = ggplot2::element_text(size = ggplot2::rel(1.14), angle = 0, vjust = 0.5),
                   axis.ticks.y = ggplot2::element_line(linewidth = 1.25, color = "black"),
                   axis.ticks.length.y = ggplot2::unit(1.25, units = "mm"),
                   axis.title.x = ggplot2::element_blank(),
                   plot.margin = ggplot2::unit(c(0.1, 0.1, 0.1, 0.1), "mm"))

  # Add x-axis style
  x_axis_theme <- switch(x_axis,
    "full" = ggplot2::theme(axis.line.x = ggplot2::element_line(linewidth = 0.3),
                            axis.ticks.x = ggplot2::element_line(color = "black"),
                            axis.ticks.length.x = ggplot2::unit(1.25, units = "mm")),
    "ticks" = ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                             axis.line.x = ggplot2::element_blank(),
                             axis.ticks.x = ggplot2::element_line(color = "black"),
                             axis.ticks.length.x = ggplot2::unit(1.25, units = "mm")),
    "blank" = ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                             axis.line.x = ggplot2::element_blank(),
                             axis.ticks.x = ggplot2::element_blank())
  )

  track_theme + x_axis_theme
}
