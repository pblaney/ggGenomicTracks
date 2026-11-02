#' @name geom_epilogo_col
#' @title Epilogos chromatin state track
#'
#' @description
#' Draws epilogos per-state scores as stacked bars spanning each bin, one color per
#' chromatin state. Positive and negative scores, as in epilogos group comparisons,
#' are stacked separately above and below zero. Data is in long format with one row
#' per bin and state; use [epilogos_long()] to convert epilogos output. `xmin`, `xmax`,
#' `y`, and `fill` default to the `start`, `end`, `score`, and `state` columns.
#'
#' Pair with [scale_fill_chromhmm()] for state colors. Also adds an "Epilogos" y-axis
#' title; add your own `labs()` after it to override it.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()].
#' @param data GRanges object or data frame of epilogos scores in long format, e.g.
#'   from [epilogos_long()].
#' @param background Fill color of the panel behind the bars, default: `"black"`. Use
#'   `NA` for none.
#' @param ... Other arguments passed to [ggplot2::layer()], e.g. `alpha`.
#' @param na.rm,show.legend,inherit.aes See [ggplot2::layer()].
#'
#' @returns List of ggplot2 layers and a y-axis title, added to a plot with `+`
#' @export
#' @keywords geom
#'
#' @examples
#' \dontrun{
#' scores <- data.table::fread("scores_mm_epilogos.chr3.txt.gz")
#'
#' ggplot2::ggplot(epilogos_long(scores)) +
#'   geom_epilogo_col() +
#'   scale_fill_chromhmm() +
#'   coord_locus("chr3:187734001-187751500", flank = 5000) +
#'   theme_genomic()
#' }
geom_epilogo_col <- function(mapping = NULL, data = NULL, background = "black", ...,
                             na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  # Default bin and state aesthetics, overridden by user mapping
  epilogo_mapping <- ggplot2::aes(xmin = start, xmax = end, y = score, fill = state)
  for(aes_name in names(mapping)) {
    epilogo_mapping[[aes_name]] <- mapping[[aes_name]]
  }

  # Background drawn as a layer rather than a theme setting so theme_genomic() keeps it
  background_layer <- NULL
  if(!is.na(background)) {
    background_layer <- ggplot2::annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf,
                                          fill = background)
  }

  list(background_layer,
       ggplot2::layer(geom = GeomEpilogoCol, stat = StatEpilogoCol, position = ggplot2::PositionStack,
                      data = as_genomic_df(data), mapping = epilogo_mapping,
                      show.legend = show.legend, inherit.aes = inherit.aes,
                      params = list(na.rm = na.rm, ...)),
       ggplot2::labs(y = "Epilogos"))
}

#' @rdname geom_epilogo_col
#' @format NULL
#' @usage NULL
#' @export
StatEpilogoCol <- ggplot2::ggproto("StatEpilogoCol", ggplot2::Stat,
  required_aes = c("xmin", "xmax", "y"),

  # Bars run from zero to the score so position_stack() stacks the states in each bin
  compute_layer = function(self, data, params, layout) {
    data$ymin <- 0
    data$ymax <- data$y
    data
  }
)

#' @rdname geom_epilogo_col
#' @format NULL
#' @usage NULL
#' @export
GeomEpilogoCol <- ggplot2::ggproto("GeomEpilogoCol", ggplot2::GeomRect,
  # Bins with 1-based closed coordinates (1-200, 201-400) leave a 1 bp gap between
  # bars; widen each by half a base on both sides so adjacent bins touch
  draw_panel = function(self, data, panel_params, coord, lineend = "butt", linejoin = "mitre") {
    data$xmin <- data$xmin - 0.5
    data$xmax <- data$xmax + 0.5
    ggplot2::GeomRect$draw_panel(data, panel_params, coord, lineend = lineend, linejoin = linejoin)
  }
)

#' @name epilogos_long
#' @title Convert epilogos scores to long format
#'
#' @description
#' Converts epilogos score output, with one score column per chromatin state, to one
#' row per bin and state with `seqnames`, `start`, `end`, `state`, and `score` columns.
#' Underscores in state names are replaced with spaces and states are ordered as their
#' columns.
#'
#' @param data Data frame or GRanges of epilogos scores. The first three columns of a
#'   data frame are taken as the chromosome, start, and end; every other column except
#'   `strand` and `width` is a state score.
#' @param state_names Names of the states in the order of the score columns, default:
#'   the 12-state ChromHMM model in epilogos order. Use `NULL` to keep the existing
#'   column names.
#'
#' @returns data.frame in long format
#' @export
#' @keywords data
epilogos_long <- function(data, state_names = chromhmm_12_states) {
  data <- as.data.frame(as_genomic_df(data))
  names(data)[1:3] <- c("seqnames", "start", "end")
  score_columns <- setdiff(names(data), c("seqnames", "start", "end", "strand", "width"))

  if(!is.null(state_names)) {
    if(length(state_names) != length(score_columns)) {
      cli::cli_abort(c("{.arg state_names} has {length(state_names)} name{?s} but {.arg data} has {length(score_columns)} score column{?s}.",
                       "i" = "Pass one name per state column, or {.code state_names = NULL} to keep the column names."))
    }
    names(data)[match(score_columns, names(data))] <- state_names
    score_columns <- state_names
  }

  long_data <- data.table::melt(data.table::as.data.table(data[, c("seqnames", "start", "end", score_columns)]),
                                id.vars = c("seqnames", "start", "end"),
                                variable.name = "state", value.name = "score")
  state_labels <- stringr::str_replace_all(score_columns, "_", " ")
  long_data$state <- factor(stringr::str_replace_all(as.character(long_data$state), "_", " "),
                            levels = state_labels)
  as.data.frame(long_data)
}

#' @name scale_fill_chromhmm
#' @title Chromatin state fill scale
#'
#' @description
#' Fill colors for the 12-state ChromHMM model used by [geom_epilogo_col()]. States are
#' matched by name with spaces, as made by [epilogos_long()], and listed in the legend
#' from promoters to heterochromatin.
#'
#' @param palette `"colorblind"` (default) for a color-blind friendly palette, or
#'   `"epilogos"` for the standard epilogos state colors.
#' @param name Legend title, default: "ChromHMM\\nStates"
#' @param ... Other arguments passed to [ggplot2::scale_fill_manual()].
#'
#' @returns ggplot2 scale
#' @export
#' @keywords scale
scale_fill_chromhmm <- function(palette = c("colorblind", "epilogos"), name = "ChromHMM\nStates", ...) {
  palette <- match.arg(palette)
  state_colors <- chromhmm_12_colors[[palette]]
  ggplot2::scale_fill_manual(name = name, values = state_colors, breaks = names(state_colors), ...)
}

# 12-state ChromHMM model, in the column order of epilogos output
chromhmm_12_states <- c("Active_Promoter", "Strong_Enhancer1", "Weak_Promoter", "Poised_Promoter",
                        "Weak_Enhancer", "Strong_Enhancer2", "Transcription_Transition",
                        "Weak_Transcription", "Transcription_Elongation", "H3K9me3_Repressed",
                        "Heterochromatin_Low_Signal", "H3K27me3_Repressed")

# State colors, in legend order from promoters to heterochromatin
chromhmm_12_colors <- list(
  colorblind = c("Active Promoter" = "#ff6db6", "Weak Promoter" = "#ffb6db", "Poised Promoter" = "#b66dff",
                 "Strong Enhancer1" = "#920000", "Strong Enhancer2" = "#924900", "Weak Enhancer" = "#db6d00",
                 "Transcription Transition" = "#006ddb", "Transcription Elongation" = "#6db6ff",
                 "Weak Transcription" = "#b6dbff", "H3K9me3 Repressed" = "#009292",
                 "H3K27me3 Repressed" = "#004949", "Heterochromatin Low Signal" = "#000000"),
  epilogos = c("Active Promoter" = "#ff0000", "Weak Promoter" = "#ffc8e6", "Poised Promoter" = "#6e1e8c",
               "Strong Enhancer1" = "#ffdc64", "Strong Enhancer2" = "#ff8200", "Weak Enhancer" = "#ffff00",
               "Transcription Transition" = "#008c64", "Transcription Elongation" = "#00ff00",
               "Weak Transcription" = "#c8ffc8", "H3K9me3 Repressed" = "#787878",
               "H3K27me3 Repressed" = "#aaaaaa", "Heterochromatin Low Signal" = "#f0f0f0")
)
