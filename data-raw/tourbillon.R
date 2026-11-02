
library(devgru)
#library(dplyr)
library(fs)
library(usethis)

# # TODO: FUNCTIONALIZE FOR ggGenomicsTracks, parametrize the number of chroms?
# # Need tick mark breaks at the midpoint of the chromosome
# linearized_hg38_end_breaks <- c()
# linearized_hg38_strip_fill <- c()
# linearized_hg38_midpoint_breaks <- c()
# for(i in 1:24) {
#   # Get end-to-end breaks for each chromosome
#   linearized_hg38_end_breaks[i] <- sum(gUtils::hg_seqlengths()[1:i])
  
#   # Set alternating strip colors for fill
#   if(i %% 2 == 1) {
#     linearized_hg38_strip_fill[i] <- "gray40"
#   } else {
#     linearized_hg38_strip_fill[i] <- "gray80"
#   }
  
#   # Get midpoint of range
#   chrom_midpoint <- round(gUtils::hg_seqlengths()[i] / 2, digits = 0)
#   # Now add the midpoints to each linear chrom after the first
#   if(i == 1) {
#     linearized_hg38_midpoint_breaks[i] <- chrom_midpoint
#   } else {
#     linearized_hg38_midpoint_breaks[i] <- sum(gUtils::hg_seqlengths()[1:i-1], chrom_midpoint)
#   }
# }
# linearized_hg38_limits <- c(1, dplyr::last(linearized_hg38_end_breaks))



# Cytoband data per chromosome
# Run from package root; reads the source file directly so the package need not be installed
cytobands_hg38 <- devgru::read_bed_file(bed_file = fs::path("inst", "extdata", "cytoband_autosome_sex_chroms.hg38.bed"),
                                        has_header = F,
                                        additional_col_names = c("name","gieStain"))
# Export for usage
usethis::use_data(cytobands_hg38, internal = FALSE, overwrite = TRUE)
