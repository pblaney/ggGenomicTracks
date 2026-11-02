#' @name label_drop_first
#' @title Blank the lowest axis label
#'
#' @description
#' Wraps a label function so the lowest visible break is left unlabeled while its tick
#' is kept. Useful in tightly stacked tracks, where the lowest y-axis label of one
#' track collides with the highest label of the track below it.
#'
#' @param labels Label function to wrap, default: [scales::label_number()]
#'
#' @returns Label function for the `labels` argument of a scale
#' @export
#' @keywords scale
#'
#' @examples
#' \dontrun{
#' ggplot2::ggplot(gistic) +
#'   geom_ridge_track(ggplot2::aes(y = G.score)) +
#'   ggplot2::scale_y_continuous(labels = label_drop_first()) +
#'   coord_locus("chr13:48303751-48481890", flank = 1e6) +
#'   theme_genomic(x_axis = "ticks")
#' }
label_drop_first <- function(labels = scales::label_number()) {
  force(labels)
  function(breaks) {
    break_labels <- labels(breaks)
    lowest <- which.min(breaks)
    if(length(lowest) > 0) {
      break_labels[lowest] <- ""
    }
    break_labels
  }
}
