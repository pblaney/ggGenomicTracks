# }}}}------->>>  Caliber ggGT-0.0.0.9000  <<<-------{{{{ #

#
# Movement Components
#

# Package-level imports live in R/ggGenomicTracks-package.R

# Set up the global default genome and number display
#' @keywords internal
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.ggGenomicTracks <- list(
    ggGenomicTracks.ref_genome = "hg38"
  )
  toset <- !(names(op.ggGenomicTracks) %in% names(op))
  if (any(toset)) options(op.ggGenomicTracks[toset])

  Sys.setenv(DEFAULT_BSGENOME = 'BSgenome.Hsapiens.UCSC.hg38::Hsapiens')
  Sys.setenv(DEFAULT_GENOME = 'BSgenome.Hsapiens.UCSC.hg38::Hsapiens')
  options(scipen = 999)

  invisible()
}

#' @name require_namespaces
#' @title Vectorised version of requireNamespace
#'
#' @description
#' Carbon copy of the helper function described by Raphael Sonabend-Friend in the package
#' `survivalmodels`. Vectorises the `requireNamespace` function and returns `TRUE` if all
#' packages, `x`, are available and `FALSE` otherwise. For internal use only.
#' https://github.com/RaphaelS1/survivalmodels/blob/main/R/helpers.R
#'
#' @param pkgs Vector of character strings with required package name(s)
#' @keywords internal
require_namespaces <- function(pkgs) {
  all(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))
}

#' @name get_data
#' @title Pulls in package data for use in package functions
#'
#' @description
#' Carbon copy of the the helper function described by stackoverflow user henfiber.
#' Pulls data into a temporary new environment and accesses it for assignment to
#' a variable in a function. For internal use only.
#' https://stackoverflow.com/questions/30951204/
#'
#' @param name_of_data Character string with name of data as it appears in package
#' @keywords internal
get_data <- function(name_of_data) {
  e <- new.env()
  name <- utils::data(list = name_of_data, package = "ggGenomicTracks", envir = e)[1]
  e[[name]]
}



##########
# DATA
##########

#' Cytoband annotations per chromosome on hg38
#'
#' Annotation of cytoband name and G-banding per chromosome
#'
#' @format A `GRanges` object with 862 ranges across chr1-22, X, and Y and 2 metadata columns:
#' \describe{
#'   \item{name}{Cytoband name, e.g. `p36.33`}
#'   \item{gieStain}{Giemsa stain result, e.g. `gneg`, `gpos25`, `acen`}
#' }
#' @source UCSC Genome Browser hg38 `cytoBand` table
#' @keywords data
"cytobands_hg38"



##########
# HELPERS
##########


# gr_absolute_ranges <- function(input_gr, cpus = 1, ref_seqlengths = gUtils::hg_seqlengths()[1:24]) {
#
#   # Set parallel cores parameter
#   message("Setting parallel cores to ", cpus, " ...")
#   doParallel::registerDoParallel(cores = cpus)
#
#   # GenomicRanges::absoluteRanges() function does not scale to 'large genomes, aka human genome
#   # Let's pararllelize the function across the chromosomes in the input GR
#   chrom_iter_list <- as.character(unique(levels(GenomeInfoDb::seqnames(input_gr))))
#
#   final_abs_ranges <- foreach::foreach(x = 1:length(chrom_iter_list), .combine = rrbind, .packages = "gUtils") %dopar% {
#
#     # Grab the ranges per chromosome
#     input_dt <- gUtils::gr2dt(input_gr)
#     ranges_per_col <- input_dt %>% dplyr::filter(seqnames == chrom_iter_list[x])
#
#     # Set base linear value depending on the chromosome
#     base_linear_val_idx <- which(chrom_iter_list[x] == names(ref_seqlengths))
#     base_linear_val <- sum(as.numeric(ref_seqlengths[1:base_linear_val_idx-1]))
#
#     # Now add the genomic coordinates to the base linear value to get the linear start and end
#     linear_start <- as.numeric(base_linear_val + ranges_per_col$start)
#     linear_end <- as.numeric(base_linear_val + ranges_per_col$end)
#
#     # Return output of the foreach loop, data.table of absolute linear coordinates
#     abs_linear_coords <- data.table::data.table("start_abs_linear" = linear_start,
#                                                 "end_abs_linear" = linear_end)
#     abs_linear_coords
#   }
#
#   # Add the absolute linear coordinates to the input GR
#   gr <- input_gr
#   gr$start_abs_linear <- final_abs_ranges$start_abs_linear
#   gr$end_abs_linear <- final_abs_ranges$end_abs_linear
#
#   # Output input GR with absolute linear coordinates
#   return(gr)
# }
#

#' @name get_linear_ranges_per_chromosome
#' @title Linearization of GRanges object for single chromosome
#'
#' @description
#' Run linearization of GRanges object where the input location relative to the individual
#' chromosome that is transformed to an absolute location on a linear end-to-end genome.
#' Useful for graphing across the whole genome.
#'
#' Note, this function was designed to be run as part of the `gr_absolute_ranges()`
#' workflow. Also, see `chrompar()` for executing in parallel.
#'
#' @param chrom_for_linearizing The chromosome to run linearization
#' @param input_ranges GenomicRanges object of interest
#' @param seq_lengths Named vector object used as template for new seq details.
#'
#' @returns data.table with absolute start and end in linear genome
#' @export
#' @keywords helper, workflow
get_linear_ranges_per_chromosome <- function(chrom_for_linearizing, input_ranges,
                                             seq_lengths) {
  # Function CLI
  cli::cli_text("{cli::symbol$pointer} {.emph {cli::col_green(chrom_for_linearizing)}}")

  # Grab the ranges per chromosome
  input_dt <- gUtils::gr2dt(input_ranges)
  ranges_per_col <- input_dt %>% dplyr::filter(seqnames == chrom_for_linearizing)

  # Set base linear value depending on the chromosome
  base_linear_val_idx <- which(chrom_for_linearizing == names(seq_lengths))
  base_linear_val <- sum(as.numeric(seq_lengths[1:base_linear_val_idx-1]))

  # Now add the genomic coordinates to the base linear value to get the linear start and end
  linear_start <- as.numeric(base_linear_val + ranges_per_col$start)
  linear_end <- as.numeric(base_linear_val + ranges_per_col$end)

  # Return data.table of absolute linear coordinates
  abs_linear_coords_per_chrom <- data.table::data.table("start_abs_linear" = linear_start,
                                                        "end_abs_linear" = linear_end)
  # Return the final imputed DT
  return(abs_linear_coords_per_chrom)
}

#' @name gr_absolute_ranges
#' @title Linearization of GRanges object
#'
#' @description
#' Run linearization of GRanges object where the input location relative to the individual
#' chromosome that is transformed to an absolute location on a linear end-to-end genome.
#' Adds 2 columns to input GRanges object:
#'   * start_abs_linear
#'   * end_abs_linear
#'
#' @param input_gr GenomicRanges object to transform
#' @param cpus Number of CPUs to use for parallel execution, default: 1
#' @param seq_lengths Named vector object used as template for new seq details.
#'
#' @returns Input GRanges object with absolute start and end in linear genome added
#' @export
#' @keywords helper, workflow
gr_absolute_ranges <- function(input_gr, cpus = 1, seq_lengths = gUtils::hg_seqlengths()[1:24]) {
  # GenomicRanges::absoluteRanges() function does not scale to large genomes, aka human genome
  # Let's pararllelize the function across the chromosomes in the input GR
  chrom_iter_list <- as.character(unique(levels(GenomeInfoDb::seqnames(input_gr))))

  final_abs_ranges <- chrompar(par_function = get_linear_ranges_per_chromosome, par_chromosomes = chrom_iter_list,
                               par_cpus = cpus, par_packages = c("gUtils", "dplyr"),
                               input_gr, seq_lengths)

  # Add the absolute linear coordinates to the input GR
  gr <- input_gr
  gr$start_abs_linear <- final_abs_ranges$start_abs_linear
  gr$end_abs_linear <- final_abs_ranges$end_abs_linear

  # Output input GR with absolute linear coordinates
  return(gr)
}


# Build the genomic coordinate x-axis for the locus of interest
# DT for continuous x-axis name, labels, limits, and breaks
get_locus_axis <- function(input_window, input_flank) {
  # Find width of desired window (native width plus flank applied to both ends)
  window_width <- GenomicRanges::width(gUtils::parse.gr(input_window) + (input_flank * 2))

  # Get labels and name
  # if total window width is >10kbp, will use Mb
  # if total window is <= 10kbp will use bp
  if(window_width > 10000) {
    window_axis_name <- paste0(stringr::str_split(string = input_window, pattern = ":")[[1]][1] %>%
                                 stringr::str_replace(pattern = "chr", replacement = "Chr "),
                               " (Mbp)")
    window_axis_labels <- scales::label_comma(accuracy = 0.001,scale = 1e-6)

  } else if (window_width <= 10000) {
    window_axis_name <- paste0(stringr::str_split(string = input_window, pattern = ":")[[1]][1] %>%
                                 stringr::str_replace(pattern = "chr", replacement = "Chr "),
                               " (bp)")
    window_axis_labels <- scales::label_comma(accuracy = 1)
  }

  # Get limits
  window_limits <- stringr::str_remove(string = input_window, pattern = "chr[\\dXY]+:") %>%
                     stringr::str_split(pattern = "-", simplify = TRUE)
  window_axis_limits <- c(as.numeric(window_limits[1,1]) - input_flank,
                          as.numeric(window_limits[1,2]) + input_flank)

  # Build the continuous genomic coordinate x-axis for the locus of interest
  locus_x_axis_params <- list("name" = window_axis_name,
                              "labels" = window_axis_labels,
                              "limits" = window_axis_limits)
  return(locus_x_axis_params)
}

locus_viewer_string <- function(input_window, input_flank = 0) {
  # Collect the parts of the viewer style string
  viewer_strings <- stringr::str_split(string = gUtils::gr.string(devgru::gr_refactor_seqs(gUtils::parse.gr(input_window) + input_flank)), pattern = ":") %>%
                      BiocGenerics::unlist() %>%
                      stringr::str_split(pattern = "-") %>%
                      BiocGenerics::unlist()

  pretty_locus_string <- paste0(viewer_strings[1],":",
                                humanize::number_as_comma(as.numeric(viewer_strings[2])),
                                "-",
                                humanize::number_as_comma(as.numeric(viewer_strings[3])),
                                "\n(",
                                humanize::number_as_comma(as.numeric(viewer_strings[3])-as.numeric(viewer_strings[2])),
                                " bp)")
  return(pretty_locus_string)
}



# intron_connetor <- function(gene_exon_dt, goi) {

#   # Find beginning and end of each gene in set
#   intron_segments <- data.table::data.table()
#   for(i in 1:length(goi)) {
#     intron_start <- gene_exon_dt[gene_name == goi[i], start][1]
#     intron_end <- gene_exon_dt[gene_name == goi[i], start][length(gene_exon_dt[gene_name == goi[i], start])]

#     # add to segment file
#     intron_segments <- gUtils::rrbind(intron_segments,
#                                       data.table::data.table("gene" = goi[i],
#                                                              "start" = intron_start,
#                                                              "end" = intron_end))
#   }
#   return(intron_segments)
# }








###########
# GEOMS
##########

#' @name geom_single_chrom_ideogrom
#' @title Single chromosome ideogram track
#'
#' @description
#' Draws the cytoband ideogram of the chromosome containing the window of interest,
#' titled with the chromosome and cytoband span, and highlights the window with a
#' zoom trace. Designed as the top track with the locus view stacked below.
#'
#' @inheritParams geom_gene_track
#' @param data Not yet implemented; the packaged `cytobands_hg38` data is used.
#' @param cytoband_file Not yet implemented; reserved for reading cytobands from file.
#' @param zoom_trace_fill_color Fill color of the window highlight and zoom trace.
#'
#' @returns ggplot object
#' @export
#' @keywords track
geom_single_chrom_ideogrom <- function(data = NULL, cytoband_file = NULL, window_of_interest = NULL,
                                       flank = 0, zoom_trace_fill_color = "#FFA319") {

  # TODO: this requies loading from scratch each time, too slow
  #       for now, let user load GTF as GRanges obj anyway
  #if(!is.null(data)) {
  #  # Quickset the data like a regular ggplot geom
  #  genes_gtf <- data
  #} else if(is.null(data) & !is.null(gtf_file)) {
  #  # Read in the GTF data
  #  genes_gtf <- devgru::read_gtf_file(gtf_file_path = gtf_file)
  #}

  # Get cytoband data
  hg38_ref_cytobands <- get_data(name_of_data = "cytobands_hg38")

  # First subset to the chromosome of WoI
  coi <- stringr::str_split(string = window_of_interest, pattern = ":")[[1]][1]
  coi_query <- hg38_ref_cytobands %Q% (seqnames == coi)

  # Generate name of view Chr CHROM (STARTCYTO - ENDCYTO)
  locus_x_axis <- get_locus_axis(input_window = window_of_interest, input_flank = flank)
  cytoband_overlap <- gUtils::gr.findoverlaps(query = devgru::gr_refactor_seqs(gUtils::parse.gr(window_of_interest) + flank),
                                              subject = hg38_ref_cytobands,
                                              scol = "name")
  # Account for small regions i.e. only 1 name
  # and for large regions, i.e. many names per cytoband covered
  if(length(cytoband_overlap) == 1) {
    cytoband_text <- cytoband_overlap$name[which.min(GenomicRanges::start(cytoband_overlap))]

  } else if(length(cytoband_overlap) > 1) {
    cytoband_text <- paste0(cytoband_overlap$name[which.min(GenomicRanges::start(cytoband_overlap))],
                            "-",
                            cytoband_overlap$name[which.max(GenomicRanges::end(cytoband_overlap))])
  }
  # Combine chromosome name and cytoband text
  ideogram_title <- stringr::str_replace(string = locus_x_axis$name, pattern = "\\(.*\\)", replacement = paste0("(", cytoband_text, ")"))

  # Build the base single chromosome ideogram
  ideogram_base <- ggplot2::ggplot() +
    ggbio::layout_karyogram(data = coi_query, cytobands = T, geom = NULL) +
    ggplot2::scale_x_continuous(position = "top", labels = scales::label_comma(accuracy = 1,scale = 1e-6, suffix = " Mb"), expand = c(0.01,0.01)) +
    ggplot2::labs(title = ideogram_title) +
    ggplot2::theme(legend.position = "none",
                   strip.background = element_blank(),
                   strip.text.y = element_blank(),
                   plot.title = element_text(size = 12, hjust = 0.5),
                   axis.line = element_blank(),
                   axis.text.x = element_text(size = 10),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_line(color = "black"),
                   axis.ticks.length.x = unit(2, units = "mm"),
                   plot.margin = unit(c(3.5,3.5,0.05,3.5), "mm"))

  # Add locus highlight and zoom trace to WoI plus flank on both ends of genomic interval of WoI
  locus_zoom_coordinates <- c(locus_x_axis$limits[1], 0, max(GenomicRanges::end(coi_query)), locus_x_axis$limits[2])

  gg_ideogram_track <- ideogram_base +
    ggplot2::geom_polygon(ggplot2::aes(x = locus_zoom_coordinates,
                                       y = c(-2, -10, -10, -2)),
                          fill = ggplot2::alpha(colour = zoom_trace_fill_color, alpha = 0.5)) +
    ggplot2::geom_rect(ggplot2::aes(xmin = locus_x_axis$limits[1],
                                    xmax = locus_x_axis$limits[2],
                                    ymin = -2,
                                    ymax = 10),
                       fill = ggplot2::alpha(colour = zoom_trace_fill_color, alpha = 0.95))

  # Output plot
  return(gg_ideogram_track)
}

#' @name geom_locus_scale
#' @title Locus coordinate and scale bar track
#'
#' @description
#' Slim track annotated with the window coordinates and width on the left and a
#' scale bar on the right, sized according to the total window width.
#'
#' @inheritParams geom_gene_track
#'
#' @returns ggplot object
#' @export
#' @keywords track
geom_locus_scale <- function(window_of_interest = NULL, flank = 0) {
  # Get the axis of the locus view
  locus_x_axis <- get_locus_axis(input_window = window_of_interest, input_flank = flank)

  # use this to make segment plot
  # By default this will have the scale on the far right
  # and will be scaled based on width of view
  loi_query <- devgru::gr_refactor_seqs(gUtils::parse.gr(window_of_interest) + flank)

  # Find the scale and appropriate abbrevation to use from total window width
  locus_scale <- dplyr::case_when(dplyr::between(x = GenomicRanges::width(loi_query), left = 50, right = 200) ~ 10,
                                  dplyr::between(x = GenomicRanges::width(loi_query), left = 201, right = 1000) ~ 40,
                                  dplyr::between(x = GenomicRanges::width(loi_query), left = 1001, right = 5000) ~ 200,
                                  dplyr::between(x = GenomicRanges::width(loi_query), left = 5001, right = 25000) ~ 1000,
                                  dplyr::between(x = GenomicRanges::width(loi_query), left = 25001, right = 100000) ~ 5000,
                                  dplyr::between(x = GenomicRanges::width(loi_query), left = 100001, right = 500000) ~ 15000,
                                  dplyr::between(x = GenomicRanges::width(loi_query), left = 500001, right = 1000000) ~ 50000,
                                  dplyr::between(x = GenomicRanges::width(loi_query), left = 1000001, right = 5000000) ~ 250000,
                                  dplyr::between(x = GenomicRanges::width(loi_query), left = 5000001, right = 20000000) ~ 500000,
                                  dplyr::between(x = GenomicRanges::width(loi_query), left = 20000001, right = 100000000) ~ 2500000,
                                  dplyr::between(x = GenomicRanges::width(loi_query), left = 100000001, right = 200000000) ~ 5000000)

  locus_scale_abbrv <- dplyr::case_when(locus_scale %in% c(10,40,200) ~ paste0(locus_scale, " bp"),
                                        locus_scale %in% c(1000,5000,15000,50000,250000,500000) ~ paste0(locus_scale / 1000, " kb"),
                                        locus_scale %in% c(2500000,5000000) ~ paste0(locus_scale / 1000000, " Mb"))

  # Build the locus scale track
  # Base plot will be barbell plot from right-hand end of the WoI
  scale_base <- ggplot2::ggplot() +
    ggplot2::geom_segment(ggplot2::aes(x = GenomicRanges::end(loi_query) - locus_scale,
                                       xend = GenomicRanges::end(loi_query) - locus_scale,
                                       y = 1.20,
                                       yend = 1.80),
                          linewidth = 0.8) +
    ggplot2::geom_segment(ggplot2::aes(x = GenomicRanges::end(loi_query),
                                       xend = GenomicRanges::end(loi_query),
                                       y = 1.20,
                                       yend = 1.80),
                          linewidth = 0.8) +
    ggplot2::geom_segment(data = gr2dt(loi_query),
                          ggplot2::aes(x = GenomicRanges::end(loi_query) - locus_scale,
                                       xend = end,
                                       y = 1.5 ,
                                       yend = 1.5),
                          linewidth = 0.8)

  # Add text annotations
  gg_scale_track <- scale_base +
    ggplot2::geom_text(ggplot2::aes(x = GenomicRanges::end(loi_query) - locus_scale,
                                    y = 1.5,
                                    label = locus_scale_abbrv),
                       size = 3.4,
                       nudge_x = locus_scale / 2,
                       nudge_y = 0.8) +
    ggplot2::geom_text(ggplot2::aes(x = GenomicRanges::start(loi_query),
                                    y = 1.5,
                                    label = locus_viewer_string(input_window = window_of_interest,
                                                                input_flank = flank)),
                       hjust = -0.03,
                       size = 3.5) +
    ggplot2::scale_x_continuous(name = locus_x_axis$name, labels = locus_x_axis$labels, limits = locus_x_axis$limits, expand = c(0.01,0.01)) +
    ggplot2::scale_y_continuous(limits = c(0,3), expand = c(0.02,0.02)) +
    ggplot2::theme(axis.line.y = element_line(linewidth = 0.3),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.05,0.05,0.05,0.05), "mm"))

  # Output plot
  return(gg_scale_track)
}

#' @name geom_gene_track
#' @title Gene structure track
#'
#' @description
#' Draws exons, UTRs, and strand-oriented intron arrows for genes of interest within
#' the window, with gene name labels placed at the 5' or 3' end.
#'
#' @param data Not yet implemented; gene structures are taken from `devgru::gene_body_hg38`.
#' @param gtf_file Not yet implemented; reserved for reading gene structures from a GTF.
#' @param genes_of_interest Character vector of gene names to draw.
#' @param window_of_interest Genomic window as a string, e.g. `"chr3:187734001-187751500"`.
#' @param flank Number of bp to extend both ends of the window, default: 0
#' @param exon_fill_color Exon fill color, or `"random"` for a distinct color per gene.
#' @param gene_name_side Side to place gene labels, `"five_prime"` or `"three_prime"`.
#'   Either a single value or one per gene in `genes_of_interest`.
#' @param gene_label_nudge Horizontal offset of gene labels in bp.
#' @param gene_label_size Gene label text size.
#' @param x_axis_theme Style of the x-axis, one of `"full"`, `"ticks"`, or `"blank"`.
#'   Use `"full"` on the bottom track of a stack and `"ticks"`/`"blank"` above it.
#'
#' @returns ggplot object
#' @export
#' @keywords track
geom_gene_track <- function(data = NULL, gtf_file = NULL, genes_of_interest,
                            window_of_interest = NULL, flank = 0, exon_fill_color = "darkred",
                            gene_name_side = "five_prime", gene_label_nudge = 100,
                            gene_label_size = 4, x_axis_theme = "full") {

  # TODO: this requies loading from scratch each time, too slow
  #       for now, let user load GTF as GRanges obj anyway
  #if(!is.null(data)) {
  #  # Quickset the data like a regular ggplot geom
  #  genes_gtf <- data
  #} else if(is.null(data) & !is.null(gtf_file)) {
  #  # Read in the GTF data
  #  genes_gtf <- devgru::read_gtf_file(gtf_file_path = gtf_file)
  #}

  # TODO: add sanity check of genes of interest

  # First subset to gene transcript structure records for query and then subset
  # to window of interest plus flank on both ends of genomic interval of WoI
  goi_query <- devgru::gene_body_hg38 %Q% (gene_name %in% c(genes_of_interest)) %&% (gUtils::parse.gr(window_of_interest) + flank)

  # Restrict the beginning and end of the GoI window if the genes extend out of the window
  # If GoI start is less than the WoI start (i.e. out-of-bounds), rrestrict start
  if(sum(GenomicRanges::start(goi_query) < GenomicRanges::start(gUtils::parse.gr(window_of_interest) + flank)) > 0) {
    goi_query <- GenomicRanges::restrict(x = goi_query, start = GenomicRanges::start(gUtils::parse.gr(window_of_interest) + flank))
  }
  # If GoI end is greater than the WoI end (i.e. out-of-bounds), restrict end
  if(sum(GenomicRanges::end(goi_query) > GenomicRanges::end(gUtils::parse.gr(window_of_interest) + flank)) > 0) {
    goi_query <- GenomicRanges::restrict(x = goi_query, end = GenomicRanges::end(gUtils::parse.gr(window_of_interest) + flank))
  }

  # Convert GoI query GR to DT
  goi_query_dt <- gUtils::gr2dt(goi_query)

  # Recode GoI strand orientation for gggenes
  goi_query_dt$orientation <- dplyr::case_when(goi_query_dt$strand == "+" ~ 1,
                                               goi_query_dt$strand == "-" ~ 0)

  # Set the aesthetic of strand orientation arrows
  intron_arrows <- grid::arrow(type = "open",
                               angle = 45,
                               length = unit(x = 3.0, units = "mm"))

  # Dynamic set of color for gene, default same color, can provide user set and
  # also a random set equal to number of genes
  if(exon_fill_color == "random") {
    exon_fill_color <- randomcoloR::distinctColorPalette(k = length(genes_of_interest))
  }

  # Build gene track plot
  # Starting with gene base - exons and UTRs
  gene_base <- ggplot2::ggplot() +
      gggenes::geom_gene_arrow(data = goi_query_dt %>% dplyr::filter(segment == "exon"),
                               ggplot2::aes(xmin = start, xmax = end, y = gene_name, fill = gene_name, forward = orientation),
                               fill = exon_fill_color,
                               color = exon_fill_color,
                               arrowhead_height = unit(9, "mm"),
                               arrowhead_width = unit(7, "mm"),
                               arrow_body_height = unit(7, "mm")) +
      ggplot2::geom_segment(data = goi_query_dt %>% dplyr::filter(segment %in% c("three_prime_utr","five_prime_utr")),
                            ggplot2::aes(x = start, xend = end, y = gene_name , yend = gene_name),
                            linewidth = 5)

  # Introns
  # To handle multiple genes with different strand orientation, add each
  # orientation type in separate layers
  # Start by pulling down the intron segments
  intron_coords <- goi_query_dt %>% dplyr::filter(segment == "intron")

  # First check if there are single or mixed orientations
  # If single, only need to add one intron arrow layer
  if(length(unique(goi_query_dt$orientation)) == 1) {
    # Check for positive strand orientation (rightward facing)
    if(unique(goi_query_dt$orientation) == 1) {
      full_gene <- gene_base +
        ggarchery::geom_arrowsegment(data = intron_coords,
                            ggplot2::aes(x = start, xend = end, y = gene_name, yend = gene_name),
                            linewidth = 0.7,
                            arrow_positions = 0.5,
                            arrows = intron_arrows)

    # Check for negative strand orientation (leftward facing)
    } else if(unique(goi_query_dt$orientation) == 0) {
      full_gene <- gene_base +
        ggarchery::geom_arrowsegment(data = intron_coords,
                            ggplot2::aes(x = end, xend = start, y = gene_name, yend = gene_name),
                            linewidth = 0.7,
                            arrow_positions = 0.5,
                            arrows = intron_arrows)
    }

  # If mixed orientation, need to add two separate intron arrow layers
  } else if (length(unique(goi_query_dt$orientation)) == 2) {
    # Add positive strand orientation (rightward facing)
    forward_gene <- gene_base +
        ggarchery::geom_arrowsegment(data = intron_coords %>% dplyr::filter(orientation == 1),
                                     ggplot2::aes(x = start, xend = end, y = gene_name, yend = gene_name),
                                     linewidth = 0.7,
                                     arrow_positions = 0.5,
                                     arrows = intron_arrows)
    # Then add the negative strand orientation (leftward facing)
    full_gene <- forward_gene +
        ggarchery::geom_arrowsegment(data = intron_coords %>% dplyr::filter(orientation == 0),
                                     ggplot2::aes(x = end, xend = start, y = gene_name, yend = gene_name),
                                     linewidth = 0.7,
                                     arrow_positions = 0.5,
                                     arrows = intron_arrows)

  }

  # Annotate with gene names
  # 5 prime vs 3 prime end
  # To handle multiple genes with both a default location for the gene name
  # and for user-defined location, link each gene of interest within the WoI
  # with the gene name side and orientation

  # First, check if there is a single side to use or if user provided set of
  # different sides to place gene names
  # TODO: need to add sanity check that if more than 1 side provided, must match length of GoI
  if(length(gene_name_side) == 1 | length(gene_name_side) == length(genes_of_interest)) {
    # Start with creating simple DT of label and y-axis discrete location (same as gene name)
    gene_name_text <- data.table("label" = genes_of_interest,
                                 "y_axis" = genes_of_interest,
                                 "side_of_anno" = paste0(gene_name_side,"_utr")) %>%
                        dplyr::left_join(y = goi_query_dt %>%
                                              dplyr::rename("label" = gene_name, "side_of_anno" = segment),
                                         by = c("label","side_of_anno")) %>%
                        dplyr::select(label,start,end,y_axis,side_of_anno,orientation)
    annotated_gene <- full_gene

    # Add each layer depending on gene name side + orientation type combination
    # Five prime side for positive strand orientation (rightward facing)
    if("five_prime_utr" %in% gene_name_text$side_of_anno & 1 %in% gene_name_text$orientation) {
      annotated_gene <- annotated_gene +
         geom_text(data = gene_name_text %>%
                            dplyr::filter(side_of_anno == "five_prime_utr" & orientation == 1),
                   aes(x = start,
                       y = y_axis,
                       label = label),
                   fontface = "italic",
                   size = gene_label_size,
                   nudge_x = -gene_label_nudge)
    }

    # Five prime side for negative strand orientation (leftward facing)
    if("five_prime_utr" %in% gene_name_text$side_of_anno & 0 %in% gene_name_text$orientation) {
      annotated_gene <- annotated_gene +
         geom_text(data = gene_name_text %>%
                            dplyr::filter(side_of_anno == "five_prime_utr" & orientation == 0),
                   aes(x = end,
                       y = y_axis,
                       label = label),
                   fontface = "italic",
                   size = gene_label_size,
                   nudge_x = gene_label_nudge)
    }

    # Three prime side for positive strand orientation (rightward facing)
    if("three_prime_utr" %in% gene_name_text$side_of_anno & 1 %in% gene_name_text$orientation) {
      annotated_gene <- annotated_gene +
         geom_text(data = gene_name_text %>%
                            dplyr::filter(side_of_anno == "three_prime_utr" & orientation == 1),
                   aes(x = end,
                       y = y_axis,
                       label = label),
                   fontface = "italic",
                   size = gene_label_size,
                   nudge_x = gene_label_nudge)
    }

    # Three prime side for negative strand orientation (leftward facing)
    if("three_prime_utr" %in% gene_name_text$side_of_anno & 0 %in% gene_name_text$orientation) {
      annotated_gene <- annotated_gene +
         geom_text(data = gene_name_text %>%
                            dplyr::filter(side_of_anno == "three_prime_utr" & orientation == 0),
                   aes(x = start,
                       y = y_axis,
                       label = label),
                   fontface = "italic",
                   size = gene_label_size,
                   nudge_x = -gene_label_nudge)
    }
  }

  # Get the settings for the locus view x-axis
  locus_x_axis <- get_locus_axis(input_window = window_of_interest, input_flank = flank)

  # Build the track
  gg_gene_track <- annotated_gene +
    ggplot2::scale_x_continuous(name = locus_x_axis$name, labels = locus_x_axis$labels, limits = locus_x_axis$limits, expand = c(0.01,0.01)) +
    ggplot2::labs(y = "Genes") +
    ggplot2::theme(legend.position = "none",
                   panel.background = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.grid.major.y = element_blank(),
                   panel.grid.major.x = element_blank(),
                   axis.line.y = element_line(linewidth = 0.3),
                   axis.title.y = element_text(size = 12.5, angle = 0, vjust = 0.5),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   plot.margin = unit(c(0.1,5,3,5), "mm"))

  # Add x-axis theme
  if(x_axis_theme == "full") {
    gg_gene_track <- gg_gene_track +
                       theme(axis.line.x = element_line(linewidth = 0.3),
                             axis.title.x = element_blank(),
                             axis.ticks.x = element_line(color = "black"),
                             axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "ticks") {
    gg_gene_track <- gg_gene_track +
                       theme(axis.title.x = element_blank(),
                             axis.text.x = element_blank(),
                             axis.line.x = element_blank(),
                             axis.ticks.x = element_line(color = "black"),
                             axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "blank") {
    gg_gene_track <- gg_gene_track +
                       theme(axis.title.x = element_blank(),
                             axis.text.x = element_blank(),
                             axis.line.x = element_blank(),
                             axis.ticks.x = element_blank())
  }

  # Output plot
  return(gg_gene_track)
}

#' @name geom_bed_track
#' @title Genomic interval track
#'
#' @description
#' Draws BED-like intervals as bars within the window, optionally colored by a
#' grouping column and sorted by width.
#'
#' @inheritParams geom_gene_track
#' @param data GRanges object of intervals.
#' @param bed_file Not yet implemented; reserved for reading intervals from a BED file.
#' @param aes_y_column Name of the column mapped to the y-axis, e.g. `"sample"`.
#' @param aes_fill_column Name of the column used to color the bars, optional.
#' @param legend_title Title of the fill legend.
#' @param aes_fill_colors Vector of colors for `aes_fill_column`; random distinct
#'   colors are used if not provided.
#' @param tornado_sort Sort rows by interval width, largest at top. Requires a
#'   `sample` column, default: FALSE
#' @param bar_height Height of the bars in mm, default: 5
#'
#' @returns ggplot object
#' @export
#' @keywords track
geom_bed_track <- function(data = NULL, bed_file = NULL, window_of_interest = NULL,
                           flank = 0, aes_y_column = NULL, aes_fill_column = NULL,
                           legend_title = NULL, aes_fill_colors = NULL, tornado_sort = F, bar_height = 5,
                           x_axis_theme = "full") {
  # TODO: multiple forms of input, like point to diff file types and use readers (VCF, MAF, etc)

  # TODO: this requies loading from scratch each time, too slow
  #       for now, let user load GTF as GRanges obj anyway
  #if(!is.null(data)) {
  #  # Quickset the data like a regular ggplot geom
  #  bed_gr <- data
  #} else if(is.null(data) & !is.null(bed_file)) {
  #  # Read in the GTF data
  #  bed_gr <- devgru::read_bed_file(gtf_file_path = bed_file)
  #}

  # First, subset to region of interest plus flank on both ends of genomic interval of WoI
  roi_query <- data %&% (gUtils::parse.gr(window_of_interest) + flank)

  # Restrict the beginning and end of the GoI window if the genes extend out of the window
  # If GoI start is less than the WoI start (i.e. out-of-bounds), rrestrict start
  if(sum(GenomicRanges::start(roi_query) < GenomicRanges::start(gUtils::parse.gr(window_of_interest) + flank)) > 0) {
    roi_query <- GenomicRanges::restrict(x = roi_query, start = GenomicRanges::start(gUtils::parse.gr(window_of_interest) + flank))
  }
  # If GoI end is greater than the WoI end (i.e. out-of-bounds), restrict end
  if(sum(GenomicRanges::end(roi_query) > GenomicRanges::end(gUtils::parse.gr(window_of_interest) + flank)) > 0) {
    roi_query <- GenomicRanges::restrict(x = roi_query, end = GenomicRanges::end(gUtils::parse.gr(window_of_interest) + flank))
  }

  # Convert RoI query GR to DT
  roi_query_dt <- gUtils::gr2dt(roi_query)

  # If there are many regions around specific locus, sort them by size with
  # largest at top of track and smallest at bottom, like a tornado
  if(tornado_sort) {
    tornado_sort_dt <- roi_query_dt %>%
                         dplyr::arrange(width)
    tornado_sort_dt$sample <- factor(tornado_sort_dt$sample,
                                     levels = unique(tornado_sort_dt$sample))
    plot_dt <- tornado_sort_dt
  } else {
    plot_dt <- roi_query_dt
  }

  # Get the settings for the locus view x-axis
  locus_x_axis <- get_locus_axis(input_window = window_of_interest, input_flank = flank)

  # Set the colors if not provided or no fill column
  if(!is.null(aes_fill_column) & !is.null(aes_fill_colors)) {
    region_fill_colors <- aes_fill_colors
  } else if(!is.null(aes_fill_column) & is.null(aes_fill_colors)) {
    # If colors not specified, set them here
    region_fill_colors <- randomcoloR::distinctColorPalette(k = plot_dt %>%
                                                                 dplyr::select(eval(aes_fill_column)) %>%
                                                                 dplyr::n_distinct())
  }

  # If provided, color bars by some grouping variable
  if(!is.null(aes_fill_column)) {
    gg_region_track <- ggplot2::ggplot(data = plot_dt,
                                       aes(xmin = start, xmax = end, y = .data[[aes_y_column]], fill = .data[[aes_fill_column]], color = .data[[aes_fill_column]])) +
                       gggenes::geom_gene_arrow(arrow_body_height = unit(bar_height, "mm"),
                                                arrowhead_height = unit(bar_height, "mm"),
                                                arrowhead_width = unit(0, "mm")) +
                       ggplot2::scale_fill_manual(name = eval(legend_title), values = region_fill_colors) +
                       ggplot2::scale_colour_manual(guide = "none", values = region_fill_colors) +
                       ggplot2::scale_x_continuous(name = locus_x_axis$name, labels = locus_x_axis$labels, limits = locus_x_axis$limits, expand = c(0.01,0.01)) +
                       ggplot2::theme(panel.background = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.grid.major.y = element_blank(),
                                      axis.line.y = element_line(linewidth = 0.3),
                                      axis.title.y = element_text(size = 12.5, angle = 0, vjust = 0.5),
                                      axis.ticks.y = element_line(linewidth = 1.25, color = "black"),
                                      axis.ticks.length.y = unit(1.25, units = "mm"),
                                      plot.margin = unit(c(0.1,0.1,0.1,0.1), "mm"))

  } else if(is.null(aes_fill_column)) {
    gg_region_track <- ggplot2::ggplot(data = plot_dt,
                                       aes(xmin = start, xmax = end, y = .data[[aes_y_column]])) +
                       gggenes::geom_gene_arrow(arrow_body_height = unit(bar_height, "mm"),
                                                arrowhead_height = unit(bar_height, "mm"),
                                                arrowhead_width = unit(0, "mm")) +
                       ggplot2::scale_x_continuous(name = locus_x_axis$name, labels = locus_x_axis$labels, limits = locus_x_axis$limits, expand = c(0.01,0.01)) +
                       ggplot2::theme(panel.background = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      panel.grid.major.y = element_blank(),
                                      axis.line.y = element_line(linewidth = 0.3),
                                      axis.title.y = element_text(size = 12.5, angle = 0, vjust = 0.5),
                                      axis.ticks.y = element_line(linewidth = 1.25, color = "black"),
                                      axis.ticks.length.y = unit(1.25, units = "mm"),
                                      plot.margin = unit(c(0.1,0.1,0.1,0.1), "mm"))
  }

  # Add x-axis theme
  if(x_axis_theme == "full") {
    gg_region_track <- gg_region_track +
                         theme(axis.line.x = element_line(linewidth = 0.3),
                               axis.title.x = element_blank(),
                               axis.ticks.x = element_line(color = "black"),
                               axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "ticks") {
    gg_region_track <- gg_region_track +
                         theme(axis.title.x = element_blank(),
                               axis.text.x = element_blank(),
                               axis.line.x = element_blank(),
                               axis.ticks.x = element_line(color = "black"),
                               axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "blank") {
    gg_region_track <- gg_region_track +
                         theme(axis.title.x = element_blank(),
                               axis.text.x = element_blank(),
                               axis.line.x = element_blank(),
                               axis.ticks.x = element_blank())
  }

  # Output plot
  return(gg_region_track)
}

# TODO: incorporate this to plotting tracks
# Remove first tick mark from y axis
# https://stackoverflow.com/questions/66279940/remove-first-label-of-y-axis
breaks_minus_one <- function(x) {
  x <- scales::pretty_breaks(n = 3)
  x[x<=0] <- NA
  x
}

#' @name geom_count_mutations
#' @title Mutation count lollipop track
#'
#' @description
#' Counts mutations at each position in the window and draws them as lollipops
#' sized by count.
#'
#' @inheritParams geom_gene_track
#' @param data data.table of mutations with `seqnames`, `start`, and `end` columns.
#' @param maf_file Not yet implemented; reserved for reading mutations from a MAF.
#' @param vcf_file Not yet implemented; reserved for reading mutations from a VCF.
#' @param aes_y_column Name of an additional column to group counts by,
#'   e.g. `"Reference_Allele"`.
#' @param aes_fill_color Fill color of the points.
#'
#' @returns ggplot object
#' @export
#' @keywords track
geom_count_mutations <- function(data = NULL, maf_file = NULL, vcf_file = NULL,
                                 window_of_interest = NULL, flank = 0,
                                 aes_y_column, aes_fill_color = "gray70",
                                 x_axis_theme = "full") {
  # TODO: multiple forms of input, like point to diff file types and use readers (VCF, MAF, etc)

  # First, subset to mutations of interest plus flank on both ends of genomic interval of WoI
  moi_query <- gr_refactor_seqs(dt2gr(data)) %&% (gUtils::parse.gr(window_of_interest) + flank)

  # Convert RoI query GR to DT
  moi_query_dt <- gUtils::gr2dt(moi_query)

  # Count the number of mutations at single locus using chromosome, start site, and reference allele
  moi_query_counts <- moi_query_dt %>%
                        dplyr::group_by(seqnames, start, .data[[aes_y_column]]) %>%
                        dplyr::count()

  # Get the settings for the locus view x-axis
  locus_x_axis <- get_locus_axis(input_window = window_of_interest, input_flank = flank)

  # Create plot
  gg_mutation_count_track <- ggplot2::ggplot(data = moi_query_counts) +
                               ggplot2::geom_segment(aes(x = start, xend = start, y = 0, yend = n), color = "gray10", alpha = 0.65, linewidth = 0.5) +
                               ggplot2::geom_point(aes(x = start, y = n, size = n), shape = 21, fill = aes_fill_color, color = aes_fill_color) +
                               ggplot2::scale_size(guide = "none") +
                               ggplot2::scale_x_continuous(name = locus_x_axis$name, labels = locus_x_axis$labels, limits = locus_x_axis$limits, expand = c(0.01,0.01)) +
                               ggplot2::scale_y_continuous(labels = scales::label_number(accuracy = 1), limits = c(0,max(moi_query_counts$n) + 1), expand = c(0.03,0.03),
                                                           breaks = scales::pretty_breaks(n = 3)) +
                               ggplot2::labs(y = "Mutation\nCount") +
                               ggplot2::theme(panel.background = element_blank(),
                                              panel.grid.major.y = element_blank(),
                                              panel.grid.major.x = element_blank(),
                                              axis.line.y = element_line(linewidth = 0.3),
                                              axis.title.y = element_text(size = 12.5, angle = 0, vjust = 0.5),
                                              axis.ticks.y = element_line(linewidth = 1.25, color = "black"),
                                              axis.ticks.length.y = unit(1.25, units = "mm"),
                                              plot.margin = unit(c(0.1,0.1,0.1,0.1), "mm"))

  # Add x-axis theme
  if(x_axis_theme == "full") {
    gg_mutation_count_track <- gg_mutation_count_track +
                                 theme(axis.line.x = element_line(linewidth = 0.3),
                                       axis.title.x = element_blank(),
                                       axis.ticks.x = element_line(color = "black"),
                                       axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "ticks") {
    gg_mutation_count_track <- gg_mutation_count_track +
                                 theme(axis.title.x = element_blank(),
                                       axis.text.x = element_blank(),
                                       axis.line.x = element_blank(),
                                       axis.ticks.x = element_line(color = "black"),
                                       axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "blank") {
    gg_mutation_count_track <- gg_mutation_count_track +
                                 theme(axis.title.x = element_blank(),
                                       axis.text.x = element_blank(),
                                       axis.line.x = element_blank(),
                                       axis.ticks.x = element_blank())
  }

  # Output plot
  return(gg_mutation_count_track)
}

#' @name geom_vaf_landscape
#' @title Variant allele frequency track
#'
#' @description
#' Scatter of mutation VAF across the window, colored by local point density.
#'
#' @inheritParams geom_gene_track
#' @param data data.table of mutations with `seqnames`, `start`, `end`, and `VAF` columns.
#' @param mut_file Not yet implemented; reserved for reading mutations from file.
#' @param vcf_file Not yet implemented; reserved for reading mutations from a VCF.
#' @param point_size Point size, default: 3
#' @param contrast Use a gray panel background, default: TRUE
#'
#' @returns ggplot object
#' @export
#' @keywords track
geom_vaf_landscape <- function(data = NULL, mut_file = NULL, vcf_file = NULL,
                               window_of_interest = NULL, flank = 0, point_size = 3,
                               contrast = T, x_axis_theme = "full") {
  # TODO: multiple forms of input, like point to diff file types and use readers (VCF, MAF, etc)

  # First, subset to mutations of interest plus flank on both ends of genomic interval of WoI
  moi_query <- gr_refactor_seqs(dt2gr(data)) %&% (gUtils::parse.gr(window_of_interest) + flank)

  # Convert RoI query GR to DT
  moi_query_dt <- gUtils::gr2dt(moi_query)

  # Get the settings for the locus view x-axis
  locus_x_axis <- get_locus_axis(input_window = window_of_interest, input_flank = flank)

  # Create plot
  gg_vaf_track <- ggplot2::ggplot(moi_query_dt) +
                    ggpointdensity::stat_pointdensity(aes(x = start, y = VAF), method = "kde2d", size = point_size) +
                    ggplot2::scale_color_gradientn(name = "Density", colours = c(paletteer::paletteer_c(palette = "grDevices::Lajolla", n = 15))) +
                    ggplot2::labs(y = "VAF") +
                    ggplot2::scale_x_continuous(name = locus_x_axis$name, labels = locus_x_axis$labels, limits = locus_x_axis$limits, expand = c(0.01,0.01)) +
                    ggplot2::scale_y_continuous(limits = c(0,1), expand = c(0.01,0.01), breaks = scales::pretty_breaks(n = 3)) +
                     ggplot2::theme(panel.grid.major.x = element_blank(),
                                    panel.grid.major.y = element_blank(),
                                    axis.line.y = element_line(linewidth = 0.3),
                                    axis.title.y = element_text(size = 12.5, angle = 0, vjust = 0.5),
                                    axis.ticks.y = element_line(linewidth = 1.25, color = "black"),
                                    axis.ticks.length.y = unit(1.25, units = "mm"),
                                    plot.margin = unit(c(0.1,0.1,0.1,0.1), "mm"))

  # Add background
  if(contrast == T) {
    gg_vaf_track <- gg_vaf_track +
                      theme(panel.background = element_rect(fill = "gray50"))

  } else if(contrast == F) {
    gg_vaf_track <- gg_vaf_track +
                      theme(panel.background = element_rect(fill = NA))
  }

  # Add x-axis theme
  if(x_axis_theme == "full") {
    gg_vaf_track <- gg_vaf_track +
                      theme(axis.line.x = element_line(linewidth = 0.3),
                            axis.title.x = element_blank(),
                            axis.ticks.x = element_line(color = "black"),
                            axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "ticks") {
    gg_vaf_track <- gg_vaf_track +
                      theme(axis.title.x = element_blank(),
                            axis.text.x = element_blank(),
                            axis.line.x = element_blank(),
                            axis.ticks.x = element_line(color = "black"),
                            axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "blank") {
    gg_vaf_track <- gg_vaf_track +
                      theme(axis.title.x = element_blank(),
                            axis.text.x = element_blank(),
                            axis.line.x = element_blank(),
                            axis.ticks.x = element_blank())
  }

  # Output plot
  return(gg_vaf_track)
}

#' @name geom_epilogo_col
#' @title Epilogos ChromHMM state track
#'
#' @description
#' Stacked bars of epilogos per-state scores across the window. Currently assumes the
#' 12-state ChromHMM model, and every column other than the coordinates is treated as
#' a state score.
#'
#' @inheritParams geom_gene_track
#' @param data data.table with `seqnames`, `start`, `end`, and one score column per state.
#' @param epilogo_file Not yet implemented; reserved for reading epilogos output from file.
#' @param contrast Use a black panel background, default: TRUE
#'
#' @returns ggplot object
#' @export
#' @keywords track
geom_epilogo_col <- function(data = NULL, epilogo_file = NULL, window_of_interest = NULL,
                             flank = 0, x_axis_theme = "full", contrast = T) {
  # TODO: # For alternate route of reading files from scratch and not already as obj
  #mm_epilogos <- aggregate_these(path_to_files = "~/morganLab/mgp1000/development/nonCoding/epiLogo/results/",
  #                             pattern_to_grab = "*_mm_epilogos.*.txt.gz",
  #                             delim = "\t",
  #                            has_header = F)
  #chromhmm_metadata <- fread(file = "nonCoding/epiLogo/metadata.tsv",
  #                         sep = "\t",
  #                         header = T)
  #colnames(mm_epilogos) <- c("chrom", "start", "end", chromhmm_metadata$long_name)

  # First, subset to epilogos of interest plus flank on both ends of genomic interval of WoI
  eoi_query <- gr_refactor_seqs(dt2gr(data)) %&% (gUtils::parse.gr(window_of_interest) + flank)

  # Restrict the beginning and end of the GoI window if the genes extend out of the window
  # If GoI start is less than the WoI start (i.e. out-of-bounds), rrestrict start
  if(sum(GenomicRanges::start(eoi_query) < GenomicRanges::start(gUtils::parse.gr(window_of_interest) + flank)) > 0) {
    eoi_query <- GenomicRanges::restrict(x = eoi_query, start = GenomicRanges::start(gUtils::parse.gr(window_of_interest) + flank))
  }
  # If GoI end is greater than the WoI end (i.e. out-of-bounds), restrict end
  if(sum(GenomicRanges::end(eoi_query) > GenomicRanges::end(gUtils::parse.gr(window_of_interest) + flank)) > 0) {
    eoi_query <- GenomicRanges::restrict(x = eoi_query, end = GenomicRanges::end(gUtils::parse.gr(window_of_interest) + flank))
  }

  # Get midpoint of EoI query GR ranges then convert to DT
  eoi_query_dt <- gUtils::gr2dt(gUtils::gr.mid(eoi_query))

  # Melt the EoI DT
  eoi_query_melt <- data.table::melt(eoi_query_dt,
                                     id.vars = c("seqnames", "start", "end"),
                                     measure.vars = colnames(eoi_query_dt)[6:length(colnames(eoi_query_dt))],
                                     variable.name = "chromhmm_state",
                                     variable.factor = T,
                                     value.name = "bits",
                                     value.factor = F)

  # Get the settings for the locus view x-axis
  locus_x_axis <- get_locus_axis(input_window = window_of_interest, input_flank = flank)

  # Remove interactions that are outside the window
  eoi_query_melt_bound <- eoi_query_melt %>%
                            dplyr::filter(start > locus_x_axis$limits[1] & end < locus_x_axis$limits[2])

  # Order the ChromHMM states
  eoi_query_melt_bound$chromhmm_state <- stringr::str_replace_all(string = eoi_query_melt_bound$chromhmm_state,
                                                                  pattern = "_",
                                                                  replacement = " ")
  eoi_query_melt_bound$chromhmm_state <- factor(eoi_query_melt_bound$chromhmm_state,
                                             levels = rev(c("Heterochromatin Low Signal","H3K27me3 Repressed","H3K9me3 Repressed",
                                                        "Weak Transcription","Transcription Elongation","Transcription Transition",
                                                        "Weak Enhancer","Strong Enhancer2","Strong Enhancer1",
                                                        "Poised Promoter","Weak Promoter","Active Promoter")))

  # Color blind friendly
  chromhmm_state_colors <- c("#ff6db6","#ffb6db","#b66dff",
                             "#920000","#924900","#db6d00",
                             "#006ddb","#6db6ff","#b6dbff",
                             "#009292","#004949","#000000")

                            #c("#260F99","#6551CC","#BFB2FF",
                            #  "#0F6B99","#51A3CC","#B2E5FF",
                            #  "#6B990F","#A3CC51","#E5FFB2",
                            #  "#99540F","#CC8E51","#FFD8B2")

  # Create plot
  gg_epilogo_track <- ggplot(eoi_query_melt_bound) +
                        ggplot2::geom_col(aes(x = start, y = bits, fill = chromhmm_state), na.rm = T, width = 200) +
                        ggplot2::labs(y = "Epilogos") +
                        scale_fill_manual(name = "ChromHMM\nStates", values = chromhmm_state_colors) +
                        ggplot2::scale_y_continuous(expand = c(0.01,0.01)) +
                        ggplot2::scale_x_continuous(name = locus_x_axis$name, labels = locus_x_axis$labels, limits = locus_x_axis$limits, expand = c(0.01,0.01)) +
                        ggplot2::theme(panel.grid.major.y = element_blank(),
                                       panel.grid.major.x = element_blank(),
                                       axis.line.y = element_line(linewidth = 0.3),
                                       axis.title.y = element_text(size = 12.5, angle = 0, vjust = 0.5),
                                       axis.ticks.y = element_blank(),
                                       axis.text.y = element_blank(),
                                       plot.margin = unit(c(0.051,0.05,0.05,0.05), "mm"))

  # Add background
  if(contrast == T) {
    gg_epilogo_track <- gg_epilogo_track +
                          theme(panel.background = element_rect(fill = "black"))

  } else if(contrast == F) {
    gg_epilogo_track <- gg_epilogo_track +
                          theme(panel.background = element_rect(fill = "white"))
  }

  # Add x-axis theme
  if(x_axis_theme == "full") {
    gg_epilogo_track <- gg_epilogo_track +
                      theme(axis.line.x = element_line(linewidth = 0.3),
                            axis.title.x = element_blank(),
                            axis.ticks.x = element_line(color = "black"),
                            axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "ticks") {
    gg_epilogo_track <- gg_epilogo_track +
                      theme(axis.title.x = element_blank(),
                            axis.text.x = element_blank(),
                            axis.line.x = element_blank(),
                            axis.ticks.x = element_line(color = "black"),
                            axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "blank") {
    gg_epilogo_track <- gg_epilogo_track +
                      theme(axis.title.x = element_blank(),
                            axis.text.x = element_blank(),
                            axis.line.x = element_blank(),
                            axis.ticks.x = element_blank())
  }

  # Output plot
  return(gg_epilogo_track)
}

#' @name geom_interactions
#' @title Genomic interaction arc track
#'
#' @description
#' Draws arcs between interacting positions that fall entirely within the window,
#' with width and color scaled by interaction score.
#'
#' @inheritParams geom_gene_track
#' @param data data.table of interactions with `seqnames`, `start`, and `end` columns.
#' @param bedpe_file Not yet implemented; reserved for reading interactions from a BEDPE.
#' @param aes_start_column Name of the column with the first anchor position.
#' @param aes_end_column Name of the column with the second anchor position.
#' @param aes_thickness_column Name of the column with the interaction score.
#' @param link_curvature Curvature of the arcs, negative values arc upward, default: -0.33
#'
#' @returns ggplot object
#' @export
#' @keywords track
geom_interactions <- function(data = NULL, bedpe_file = NULL, window_of_interest = NULL,
                              flank = 0, aes_start_column, aes_end_column, aes_thickness_column,
                              link_curvature = -0.33, x_axis_theme = "full") {
  # TODO: multiple forms of input, like point to diff file types and use readers (VCF, MAF, etc)

  # First, subset to interactions of interest plus flank on both ends of genomic interval of WoI
  ioi_query <- gr_refactor_seqs(dt2gr(data)) %&% (gUtils::parse.gr(window_of_interest) + flank)

  # Convert RoI query GR to DT
  ioi_query_dt <- gUtils::gr2dt(ioi_query)

  # Create a DT of necessary columns as input for MintChip
  ioi_mintchip_dt <- ioi_query_dt %>%
                       dplyr::select(all_of(c(aes_start_column, aes_end_column, aes_thickness_column))) %>%
                       dplyr::rename("start" = all_of(aes_start_column),
                                     "end" = all_of(aes_end_column),
                                     "thickness" = all_of(aes_thickness_column))

  # Need to convert integer format to numeric
  ioi_mintchip_dt$start <- as.numeric(ioi_mintchip_dt$start)
  ioi_mintchip_dt$end <- as.numeric(ioi_mintchip_dt$end)

  # Get the settings for the locus view x-axis
  locus_x_axis <- get_locus_axis(input_window = window_of_interest, input_flank = flank)

  # Remove interactions that are outside the window
  ioi_mintchip_dt <- ioi_mintchip_dt %>%
                       dplyr::filter(start > locus_x_axis$limits[1] & end < locus_x_axis$limits[2])

  # Create plot
  gg_interaction_track <- ggplot(data = ioi_mintchip_dt) +
      ggplot2::geom_curve(aes(x = start, xend = end, y = 0, yend = 0, linewidth = thickness, color = thickness),
                           lineend = "round",
                           alpha = 0.8,
                           curvature = link_curvature,
                           ncp = 20) +
      ggplot2::scale_color_gradientn(name = "EPIraction\nScore", colours = c(paletteer::paletteer_c("grDevices::Berlin", 30)[15:30])) +
      ggplot2::scale_linewidth(guide = NULL) +
      ggplot2::scale_x_continuous(name = locus_x_axis$name, labels = locus_x_axis$labels, limits = locus_x_axis$limits, expand = c(0.01,0.01)) +
      ggplot2::scale_y_continuous(limits = c(0,0.01), expand = c(0,0)) +
      ggplot2::labs(y = "Interactions") +
      ggplot2::theme(panel.background = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.grid.major.y = element_blank(),
                     axis.line.y = element_line(linewidth = 0.3),
                     axis.title.y = element_text(size = 14),
                     axis.ticks.y = element_blank(),
                     axis.text.y = element_blank(),
                     plot.margin = unit(c(0.1,0.1,0.1,0.1), "mm"))

  # Add x-axis theme
  if(x_axis_theme == "full") {
    gg_interaction_track <- gg_interaction_track +
                              theme(axis.line.y = element_line(linewidth = 0.3),
                                    axis.title.x = element_blank(),
                                    axis.ticks.x = element_line(color = "black"),
                                    axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "ticks") {
    gg_interaction_track <- gg_interaction_track +
                              theme(axis.title.x = element_blank(),
                                    axis.text.x = element_blank(),
                                    axis.line.x = element_blank(),
                                    axis.ticks.x = element_line(color = "black"),
                                    axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "blank") {
    gg_interaction_track <- gg_interaction_track +
                              theme(axis.title.x = element_blank(),
                                    axis.text.x = element_blank(),
                                    axis.line.x = element_blank(),
                                    axis.ticks.x = element_blank())
  }

  # Output plot
  return(gg_interaction_track)
}

# add functionality to connect in-window DEL, INV, and DUP BPs
# but maintain distinct SV bps
#' @name geom_sv_breakpoints
#' @title Structural variant breakpoint track
#'
#' @description
#' Draws SV breakpoints within the window as points colored by SV class
#' (`DEL`, `DUP`, `INV`, `TRA`).
#'
#' @inheritParams geom_gene_track
#' @param data data.table of breakpoints with `seqnames`, `start`, and `end` columns.
#' @param bedpe_file Not yet implemented; reserved for reading breakpoints from a BEDPE.
#' @param aes_y_column Name of the column mapped to the y-axis.
#' @param aes_fill_column Name of the column with the SV class.
#' @param aes_shape_column Not yet implemented; reserved for mapping point shape.
#'
#' @returns ggplot object
#' @export
#' @keywords track
geom_sv_breakpoints <- function(data = NULL, bedpe_file = NULL, window_of_interest = NULL,
                                flank = 0, aes_y_column, aes_fill_column, aes_shape_column,
                                x_axis_theme = "full") {
  # TODO: # For alternate route of reading files from scratch and not already as obj
  #mm_epilogos <- aggregate_these(path_to_files = "~/morganLab/mgp1000/development/nonCoding/epiLogo/results/",
  #                             pattern_to_grab = "*_mm_epilogos.*.txt.gz",
  #                             delim = "\t",
  #                            has_header = F)
  #chromhmm_metadata <- fread(file = "nonCoding/epiLogo/metadata.tsv",
  #                         sep = "\t",
  #                         header = T)
  #colnames(mm_epilogos) <- c("chrom", "start", "end", chromhmm_metadata$long_name)

   # First, subset to breakpoints of interest plus flank on both ends of genomic interval of WoI
  boi_query <- gr_refactor_seqs(dt2gr(data)) %&% (gUtils::parse.gr(window_of_interest) + flank)

  # Convert BoI query GR to DT
  boi_query_dt <- gUtils::gr2dt(boi_query)

  # Get the settings for the locus view x-axis
  locus_x_axis <- get_locus_axis(input_window = window_of_interest, input_flank = flank)

  # Set the color palette
  sv_type_colors <- data.table::data.table("sv_type" = c("DEL","DUP","INV","TRA"),
                                           "colors" = paletteer::paletteer_d("ggsci::default_jama")[1:4])

  # TODO: for different y-axis variables need to account for when setting y-axis scale
  # also account for use of color when provided shape, while fill used for no shape
  # Create plot
  gg_breakpoint_track <- ggplot2::ggplot(data = boi_query_dt) +
                     ggplot2::geom_point(aes(x = start,
                                             y = .data[[aes_y_column]],
                                             fill = .data[[aes_fill_column]],
                                             color = .data[[aes_fill_column]]),
                                         shape = 21, size = 4) +
                     scale_fill_manual(name = "SV\nClass",
                                       values = sv_type_colors %>%
                                                  dplyr::filter(sv_type %in% unique(boi_query_dt %>%
                                                                                      dplyr::select(dplyr::all_of(aes_fill_column)) %>%
                                                                                      purrr::as_vector() %>%
                                                                                      unique())) %>%
                                                                  dplyr::select(colors) %>%
                                                                  purrr::as_vector() %>%
                                                                  as.character()) +
                    scale_color_manual(guide = "none",
                                      values = sv_type_colors %>%
                                                dplyr::filter(sv_type %in% unique(boi_query_dt %>%
                                                                                   dplyr::select(dplyr::all_of(aes_fill_column)) %>%
                                                                                   purrr::as_vector() %>%
                                                                                   unique())) %>%
                                                dplyr::select(colors) %>%
                                                purrr::as_vector() %>%
                                                as.character()) +
                     ggplot2::labs(y = "Breakpoints") +
                     ggplot2::scale_y_discrete(expand = c(0.03,0.03)) +
                     ggplot2::scale_x_continuous(name = locus_x_axis$name, labels = locus_x_axis$labels, limits = locus_x_axis$limits, expand = c(0.01,0.01)) +
                      ggplot2::theme(legend.position = "right",
                                     panel.background = element_blank(),
                                     panel.grid.major.y = element_blank(),
                                     panel.grid.major.x = element_blank(),
                                     axis.line.y = element_line(linewidth = 0.3),
                                     axis.ticks.y = element_blank(),
                                     axis.ticks.length.y = unit(1.25, units = "mm"),
                                     axis.text.y = element_blank(),
                                     axis.title.y = element_text(size = 12.5, angle = 0, vjust = 0.5),
                                     plot.margin = unit(c(0.05,0.05,0.05,0.05), "mm"))

  # Add x-axis theme
  if(x_axis_theme == "full") {
    gg_breakpoint_track <- gg_breakpoint_track +
                            theme(axis.line.x = element_line(linewidth = 0.3),
                                  axis.title.x = element_blank(),
                                  axis.ticks.x = element_line(color = "black"),
                                  axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "ticks") {
    gg_breakpoint_track <- gg_breakpoint_track +
                            theme(axis.title.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.line.x = element_blank(),
                                  axis.ticks.x = element_line(color = "black"),
                                  axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "blank") {
    gg_breakpoint_track <- gg_breakpoint_track +
                            theme(axis.title.x = element_blank(),
                                  axis.text.x = element_blank(),
                                  axis.line.x = element_blank(),
                                  axis.ticks.x = element_blank())
  }

  # Output plot
  return(gg_breakpoint_track)
}

# TODO: take bigwig/bedgraph-like/extended BED like
# designed to take files that look like:
#chr     start   end     MM59300 MM58698 ... MM54
#chr1    10052   10447   3.36224889424292        5.15650463271808
#chr1    180676  181211  3.77952639507682        5.03959938
#' @name geom_ridge_track
#' @title Signal track
#'
#' @description
#' Draws a quantitative signal (e.g. bedGraph-like coverage) as bars at the
#' midpoint of each range within the window.
#'
#' @inheritParams geom_gene_track
#' @param data GRanges object with one or more signal columns.
#' @param file Not yet implemented; reserved for reading bigWig/bedGraph files.
#' @param aes_y_column Name of the signal column, or the name given to the melted
#'   value column when `aes_facet_y_column` is used.
#' @param aes_facet_y_column Name given to the melted sample column when `data`
#'   has one signal column per sample, optional.
#' @param fill_color Fill color of the bars.
#'
#' @returns ggplot object
#' @export
#' @keywords track
geom_ridge_track <- function(data = NULL, file = NULL, window_of_interest = NULL,
                             flank = 0, aes_y_column = NULL, aes_facet_y_column = NULL,
                             fill_color = "darkred", x_axis_theme = "full") {
  # TODO: multiple forms of input, like point to diff file types and use readers (bigWig, bedGraph, etc)

  # TODO: this requies loading from scratch each time, too slow
  #       for now, let user load GTF as GRanges obj anyway
  #if(!is.null(data)) {
  #  # Quickset the data like a regular ggplot geom
  #  bed_gr <- data
  #} else if(is.null(data) & !is.null(bed_file)) {
  #  # Read in the GTF data
  #  bed_gr <- devgru::read_bed_file(gtf_file_path = bed_file)
  #}

  # First, subset to region of interest plus flank on both ends of genomic interval of WoI
  roi_query <- data %&% (gUtils::parse.gr(window_of_interest) + flank)

  # Restrict the beginning and end of the GoI window if the genes extend out of the window
  # If GoI start is less than the WoI start (i.e. out-of-bounds), rrestrict start
  if(sum(GenomicRanges::start(roi_query) < GenomicRanges::start(gUtils::parse.gr(window_of_interest) + flank)) > 0) {
    roi_query <- GenomicRanges::restrict(x = roi_query, start = GenomicRanges::start(gUtils::parse.gr(window_of_interest) + flank))
  }
  # If GoI end is greater than the WoI end (i.e. out-of-bounds), restrict end
  if(sum(GenomicRanges::end(roi_query) > GenomicRanges::end(gUtils::parse.gr(window_of_interest) + flank)) > 0) {
    roi_query <- GenomicRanges::restrict(x = roi_query, end = GenomicRanges::end(gUtils::parse.gr(window_of_interest) + flank))
  }

  # Find midpoint of range for each value
  roi_mid_query <- gUtils::gr.mid(roi_query)

  # Convert RoI query GR to DT
  roi_query_dt <- gUtils::gr2dt(roi_mid_query)

  # Get the settings for the locus view x-axis
  locus_x_axis <- get_locus_axis(input_window = window_of_interest, input_flank = flank)

  # Need to melt the RoI query DT get into ggplot format with multiple samples
  # with columns for scores
  if(!is.null(aes_facet_y_column)) {
    melt_roi_query_dt <- data.table::melt(data = roi_query_dt[,-c("strand", "end", "width")],
                                          id.vars = c("seqnames", "start"),
                                          value.name = eval(aes_y_column),
                                          variable.name = aes_facet_y_column)
  } else if(is.null(aes_facet_y_column)){
    melt_roi_query_dt <- roi_query_dt %>%
                          dplyr::select(seqnames,start,eval(aes_y_column))
  }

  gg_ridge_plot <- ggplot2::ggplot(data = melt_roi_query_dt,
                                     aes(x = start, y = .data[[aes_y_column]])) +
                      ggplot2::geom_col(fill = fill_color, width = 150) +
                      ggplot2::scale_x_continuous(name = locus_x_axis$name, labels = locus_x_axis$labels, limits = locus_x_axis$limits, expand = c(0.01,0.01)) +
                      ggplot2::scale_y_continuous(breaks = scales::pretty_breaks(n = 3), expand = c(0.02,0.02)) +
                      ggplot2::theme(panel.background = element_blank(),
                                     panel.grid.minor = element_blank(),
                                     panel.grid.major.y = element_blank(),
                                     axis.line.y = element_line(linewidth = 0.3),
                                     axis.title.y = element_text(size = 12.5, angle = 0, vjust = 0.5),
                                     axis.ticks.y = element_line(linewidth = 1.25, color = "black"),
                                     axis.ticks.length.y = unit(1.25, units = "mm"),
                                     plot.margin = unit(c(0.1,0.1,0.1,0.1), "mm"))
    # Add x-axis theme
  if(x_axis_theme == "full") {
    gg_ridge_plot <- gg_ridge_plot +
                         theme(axis.line.x = element_line(linewidth = 0.3),
                               axis.title.x = element_blank(),
                               axis.ticks.x = element_line(color = "black"),
                               axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "ticks") {
    gg_ridge_plot <- gg_ridge_plot +
                         theme(axis.title.x = element_blank(),
                               axis.text.x = element_blank(),
                               axis.line.x = element_blank(),
                               axis.ticks.x = element_line(color = "black"),
                               axis.ticks.length.x = unit(1.25, units = "mm"))

  } else if(x_axis_theme == "blank") {
    gg_ridge_plot <- gg_ridge_plot +
                         theme(axis.title.x = element_blank(),
                               axis.text.x = element_blank(),
                               axis.line.x = element_blank(),
                               axis.ticks.x = element_blank())
  }

  # first, zoom to region provided
  # TODO: can take input DT or GR
  #zoom_bed_gr <- data %&% (gUtils::parse.gr(window) + pad)
  # Find midpoint of range for each value
  #zoom_mid_bed_gr <- gUtils::gr.mid(zoom_bed_gr)
  # Convert to data.table
  #zoom_mid_bed_dt <- gUtils::gr2dt(zoom_mid_bed_gr)

  # Need to melt to get into ggplot format
  #melt_zoom_bed_dt <- data.table::melt(data = zoom_mid_bed_dt[,-c("strand", "width")],
  #                                     id.vars = c("seqnames", "start", "end"),
  #                                     value.name = aes_y_axis_column, # TODO: PARAMETERIZE
  #                                     variable.name = aes_facet_y_column) # TODO: PARAMETRIZE

  #gg_track_plot <- ggplot(data = melt_zoom_bed_dt) +
  #                   ggHoriPlot::geom_horizon(aes(x = start, y = .data[[aes_y_axis_column]], fill = ..Cutpoints..),
  #                                            origin = 0,
  #                                            horizonscale = c(0,0.01, round(max(melt_zoom_bed_dt[[aes_y_axis_column]]), digits = 0))) +
  #                   scale_fill_manual(values = c(aes_fill_color, "transparent")) +
  #                   ggplot2::scale_x_continuous(labels = scales::label_bytes(accuracy = 0.001), limits = c(get_window_limits(input_window = window)$left - pad, get_window_limits(input_window = window)$right + pad), expand = c(0.01, 0.01)) +
  #                   scale_y_continuous(n.breaks = 3,
  #                                      labels = scales::label_number(scale = round(max(melt_zoom_bed_dt[[aes_y_axis_column]]) + 1, digits = 0),
  #                                                                    accuracy = 1)) +
  #                   facet_grid(.data[[aes_facet_y_column]] ~ ., switch = "y") +
  #                   theme(legend.position = "none",
  #                         panel.background = element_rect(fill = NA),
  #                         panel.grid.major.y = element_line(color = "gray80"),
  #                         axis.title.x = element_blank(),
  #                         axis.line.x = element_blank(),
  #                         axis.text.x = element_blank(),
  #                         axis.ticks.x = element_line(color = "black"),
  #                         axis.ticks.length.x = unit(3, units = "mm"),
  #                         panel.spacing.y = unit(1, "lines"),
  #                         strip.text.y.left = element_text(size = 14, angle = 0, hjust = 0),
  #                         strip.placement = "outside",
  #                         strip.background = element_rect(fill = NA),
  #                         axis.ticks.y = element_line(linewidth = 3, color = "black"),
  #                         axis.ticks.length.y = unit(2, units = "mm"),
  #                         axis.text.y = element_text(size = 14, face = "bold"),
  #                         axis.title.y = element_blank())

  return(gg_ridge_plot)
}


