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
