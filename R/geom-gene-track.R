#' @name geom_gene_track
#' @title Gene structure track
#'
#' @description
#' Draws gene models with exons, UTRs, and introns with strand direction arrows, plus
#' one italic label per gene at its 5' or 3' end. One row per level of the `y`
#' aesthetic, which defaults to `gene_name`. Pair with [coord_locus()] to zoom to a
#' window and `theme_genomic(y_axis = "blank")` for a clean gene track.
#'
#' Labels for genes extending past the window are kept at the panel edge.
#'
#' @param mapping Set of aesthetic mappings created by [ggplot2::aes()]. Defaults to
#'   `y = gene_name`. Map `fill` (e.g. `fill = gene_name`) to color exons by a column.
#' @param data GRanges object or data frame of gene segments with `seqnames`, `start`,
#'   `end`, `strand`, `gene_name`, and `segment` (`exon`, `intron`, `five_prime_utr`,
#'   `three_prime_utr`) columns. Defaults to `devgru::gene_body_hg38`.
#' @param genes_of_interest Character vector of gene names to draw. Default `NULL`
#'   draws every gene in the window.
#' @param exon_fill_color Exon fill color, used when `fill` is not mapped.
#' @param gene_name_side Side to place gene labels, `"five_prime"` or `"three_prime"`.
#'   Either a single value or one per gene in `genes_of_interest`.
#' @param gene_label_size Gene label text size, default: 4
#' @param gene_label_padding Space between gene end and label in mm, default: 1.5
#' @param exon_height Height of exons in mm, default: 7
#' @param exon_arrowhead_width Width of exon arrowheads in mm. Set to 0 for plain
#'   rectangular exons, default: 7
#' @param utr_linewidth Line width of UTRs, default: 5
#' @param intron_linewidth Line width of introns, default: 0.7
#' @param inherit.aes See [ggplot2::layer()].
#'
#' @returns List of ggplot2 layers
#' @export
#' @keywords geom
#'
#' @examples
#' \dontrun{
#' ggplot2::ggplot() +
#'   geom_gene_track(genes_of_interest = "BCL6") +
#'   coord_locus("chr3:187734001-187751500", flank = 5000) +
#'   theme_genomic(y_axis = "blank")
#'
#' # All genes in the window, exons colored per gene
#' ggplot2::ggplot() +
#'   geom_gene_track(ggplot2::aes(fill = gene_name)) +
#'   coord_locus("chr3:187600001-187800000") +
#'   theme_genomic(y_axis = "blank")
#' }
geom_gene_track <- function(mapping = NULL, data = NULL, genes_of_interest = NULL,
                            exon_fill_color = "darkred", gene_name_side = "five_prime",
                            gene_label_size = 4, gene_label_padding = 1.5,
                            exon_height = 7, exon_arrowhead_width = 7,
                            utr_linewidth = 5, intron_linewidth = 0.7, inherit.aes = FALSE) {
  # Default gene structure data
  if(is.null(data)) {
    data <- devgru::gene_body_hg38
  }
  gene_dt <- prep_gene_segments(data = data, genes_of_interest = genes_of_interest)

  # Row aesthetic shared by every layer, overridden by user mapping
  y_aes <- if(!is.null(mapping$y)) mapping$y else ggplot2::aes(y = gene_name)$y
  user_aes <- mapping[setdiff(names(mapping), "y")]

  # Exons: user fill mapping colors both fill and outline, otherwise fixed color
  exon_mapping <- ggplot2::aes(xmin = start, xmax = end, forward = forward)
  exon_mapping$y <- y_aes
  for(aes_name in names(user_aes)) {
    exon_mapping[[aes_name]] <- user_aes[[aes_name]]
  }
  exon_params <- list()
  if(is.null(exon_mapping$fill)) {
    exon_params <- list(fill = exon_fill_color, colour = exon_fill_color)
  } else if(is.null(exon_mapping$colour)) {
    exon_mapping$colour <- exon_mapping$fill
  }
  exon_layer <- do.call(gggenes::geom_gene_arrow,
                        c(list(mapping = exon_mapping,
                               data = gene_dt[gene_dt$segment == "exon", ],
                               arrowhead_height = ggplot2::unit(exon_height + 2, "mm"),
                               arrowhead_width = ggplot2::unit(exon_arrowhead_width, "mm"),
                               arrow_body_height = ggplot2::unit(exon_height, "mm"),
                               inherit.aes = inherit.aes),
                          exon_params))

  # UTRs as thick segments
  utr_mapping <- ggplot2::aes(x = start, xend = end)
  utr_mapping$y <- y_aes
  utr_mapping$yend <- y_aes
  utr_layer <- ggplot2::geom_segment(mapping = utr_mapping,
                                     data = gene_dt[gene_dt$segment %in% c("five_prime_utr", "three_prime_utr"), ],
                                     linewidth = utr_linewidth, inherit.aes = inherit.aes)

  # Introns as lines with an arrow at the midpoint pointing in the direction of transcription
  intron_mapping <- ggplot2::aes(x = arrow_from, xend = arrow_to)
  intron_mapping$y <- y_aes
  intron_mapping$yend <- y_aes
  intron_layer <- ggarchery::geom_arrowsegment(mapping = intron_mapping,
                                               data = gene_dt[gene_dt$segment == "intron", ],
                                               linewidth = intron_linewidth,
                                               arrow_positions = 0.5,
                                               arrows = grid::arrow(type = "open", angle = 45,
                                                                    length = ggplot2::unit(3, "mm")),
                                               inherit.aes = inherit.aes)

  # One label per gene at the requested end
  label_mapping <- ggplot2::aes(x = label_x, xmin = start, xmax = end, label = gene_name, hjust = label_hjust)
  label_mapping$y <- y_aes
  label_layer <- ggplot2::layer(geom = GeomGeneLabel, stat = ggplot2::StatIdentity,
                                position = ggplot2::PositionIdentity,
                                data = get_gene_labels(gene_dt = gene_dt, genes_of_interest = genes_of_interest,
                                                       gene_name_side = gene_name_side),
                                mapping = label_mapping, inherit.aes = inherit.aes, show.legend = FALSE,
                                params = list(size = gene_label_size, fontface = "italic",
                                              padding = gene_label_padding,
                                              raise = exon_height / 2 + 1.5))

  list(intron_layer, utr_layer, exon_layer, label_layer)
}

#' @name prep_gene_segments
#' @title Prepare gene segment data for plotting
#'
#' @description
#' Converts to a data frame, subsets to genes of interest, and adds strand direction
#' columns used by the exon and intron layers. For internal use only.
#'
#' @inheritParams geom_gene_track
#' @keywords internal
prep_gene_segments <- function(data, genes_of_interest = NULL) {
  gene_dt <- as_genomic_df(data)

  if(!is.null(genes_of_interest)) {
    missing_genes <- setdiff(genes_of_interest, gene_dt$gene_name)
    if(length(missing_genes) > 0) {
      cli::cli_warn("Gene{?s} not found in gene data: {.val {missing_genes}}")
    }
    gene_dt <- gene_dt[gene_dt$gene_name %in% genes_of_interest, ]
  }

  # Strand direction: exon arrows point forward on + strand, introns drawn 5' to 3'
  is_reverse <- as.character(gene_dt$strand) == "-"
  gene_dt$forward <- !is_reverse
  gene_dt$arrow_from <- ifelse(is_reverse, gene_dt$end, gene_dt$start)
  gene_dt$arrow_to <- ifelse(is_reverse, gene_dt$start, gene_dt$end)
  gene_dt
}

#' @name get_gene_labels
#' @title One label position per gene
#'
#' @description
#' Collapses gene segments to one row per gene spanning the whole gene, with the label
#' anchored at the 5' or 3' end and justified away from the gene. For internal use only.
#'
#' @param gene_dt Data frame from `prep_gene_segments()`
#' @inheritParams geom_gene_track
#' @keywords internal
get_gene_labels <- function(gene_dt, genes_of_interest = NULL, gene_name_side = "five_prime") {
  if(nrow(gene_dt) == 0) {
    return(gene_dt)
  }

  # Keep the first row per gene for any other mapped columns, spanning the whole gene
  gene_labels <- gene_dt[!duplicated(gene_dt$gene_name), ]
  gene_labels$start <- as.vector(tapply(gene_dt$start, gene_dt$gene_name, min)[gene_labels$gene_name])
  gene_labels$end <- as.vector(tapply(gene_dt$end, gene_dt$gene_name, max)[gene_labels$gene_name])

  # Side per gene, either shared or matched to genes of interest
  if(length(gene_name_side) == 1) {
    side <- rep(gene_name_side, nrow(gene_labels))
  } else if(!is.null(genes_of_interest) && length(gene_name_side) == length(genes_of_interest)) {
    side <- gene_name_side[match(gene_labels$gene_name, genes_of_interest)]
  } else {
    cli::cli_abort("{.arg gene_name_side} must be length 1 or match the length of {.arg genes_of_interest}.")
  }
  if(!all(side %in% c("five_prime", "three_prime"))) {
    cli::cli_abort("{.arg gene_name_side} must be {.val five_prime} or {.val three_prime}.")
  }

  # Label at the left end of the gene sits left of it, and vice versa
  label_at_left <- (side == "five_prime") == gene_labels$forward
  gene_labels$label_x <- ifelse(label_at_left, gene_labels$start, gene_labels$end)
  gene_labels$label_hjust <- ifelse(label_at_left, 1, 0)
  gene_labels
}

#' @name place_gene_labels
#' @title Position gene labels within the panel at draw time
#'
#' @description
#' Places each label at its requested gene end if it fits in the panel, otherwise at
#' the opposite gene end, otherwise centered above the visible part of the gene.
#' Called by the label grob once the panel size is known. For internal use only.
#'
#' @param label_grob gTree of class `gene_label_grob` built by `GeomGeneLabel`
#' @keywords internal
place_gene_labels <- function(label_grob) {
  lab <- label_grob$label_data

  # Label widths and padding in panel npc, known only once the panel is drawn
  label_width <- vapply(seq_len(nrow(lab)), function(i) {
    grid::convertWidth(grid::grobWidth(grid::textGrob(lab$label[i], gp = grid::gpar(fontsize = lab$fontsize[i],
                                                                                    fontface = lab$fontface[i],
                                                                                    fontfamily = lab$family[i]))),
                       "npc", valueOnly = TRUE)
  }, numeric(1))
  pad <- grid::convertWidth(grid::unit(label_grob$padding, "mm"), "npc", valueOnly = TRUE)

  # A label beside a gene end fits if its anchor is visible and the text stays in the panel
  fits_beside <- function(anchor, hjust) {
    ifelse(hjust == 1,
           anchor <= 1 & anchor - pad - label_width >= 0,
           anchor >= 0 & anchor + pad + label_width <= 1)
  }
  alt_x <- ifelse(lab$hjust == 1, lab$gene_end, lab$gene_start)
  alt_hjust <- 1 - lab$hjust
  use_main <- fits_beside(lab$x, lab$hjust)
  use_alt <- !use_main & fits_beside(alt_x, alt_hjust)
  use_above <- !use_main & !use_alt

  # Above: centered on the visible part of the gene, kept within the panel
  visible_mid <- (pmax(lab$gene_start, 0) + pmin(lab$gene_end, 1)) / 2
  above_x <- pmin(pmax(visible_mid, label_width / 2), 1 - label_width / 2)

  x <- ifelse(use_main, lab$x, ifelse(use_alt, alt_x, above_x))
  hjust <- ifelse(use_main, lab$hjust, ifelse(use_alt, alt_hjust, 0.5))
  x_pad_mm <- ifelse(use_above, 0, ifelse(hjust == 1, -1, 1) * label_grob$padding)
  y_raise_mm <- ifelse(use_above, label_grob$raise, 0)

  label_text <- grid::textGrob(label = lab$label,
                               x = grid::unit(x, "npc") + grid::unit(x_pad_mm, "mm"),
                               y = grid::unit(lab$y, "npc") + grid::unit(y_raise_mm, "mm"),
                               hjust = hjust, vjust = ifelse(use_above, 0, 0.5),
                               gp = grid::gpar(col = lab$col, fontsize = lab$fontsize,
                                               fontface = lab$fontface, fontfamily = lab$family))
  grid::setChildren(label_grob, grid::gList(label_text))
}

#' @exportS3Method grid::makeContent
makeContent.gene_label_grob <- function(x) {
  place_gene_labels(x)
}

#' @rdname geom_gene_track
#' @format NULL
#' @usage NULL
#' @export
GeomGeneLabel <- ggplot2::ggproto("GeomGeneLabel", ggplot2::GeomText,
  optional_aes = c("xmin", "xmax"),

  draw_panel = function(data, panel_params, coord, padding = 1.5, raise = 5, na.rm = FALSE) {
    # Positions to panel npc, including the gene span in xmin/xmax
    data <- coord$transform(data, panel_params)

    # Final placement needs label widths in the drawn panel, see place_gene_labels()
    grid::gTree(label_data = data.frame(label = data$label, x = data$x, y = data$y, hjust = data$hjust,
                                        gene_start = data$xmin, gene_end = data$xmax,
                                        col = ggplot2::alpha(data$colour, data$alpha),
                                        fontsize = data$size * ggplot2::.pt,
                                        fontface = data$fontface, family = data$family),
                padding = padding, raise = raise, cl = "gene_label_grob")
  }
)
