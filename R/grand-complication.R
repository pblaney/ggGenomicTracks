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

#' @name legacy_track_params
#' @title Shared parameters of tracks not yet converted to layers
#'
#' @param window_of_interest Genomic window as a string, e.g. `"chr3:187734001-187751500"`.
#' @param flank Number of bp to extend both ends of the window, default: 0
#' @param x_axis_theme Style of the x-axis, one of `"full"`, `"ticks"`, or `"blank"`.
#'   Use `"full"` on the bottom track of a stack and `"ticks"`/`"blank"` above it.
#' @keywords internal
NULL


#' @name geom_interactions
#' @title Genomic interaction arc track
#'
#' @description
#' Draws arcs between interacting positions that fall entirely within the window,
#' with width and color scaled by interaction score.
#'
#' @inheritParams legacy_track_params
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
#' @inheritParams legacy_track_params
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


