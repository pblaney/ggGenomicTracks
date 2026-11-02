# ggGenomicTracks

Genomic track visualizations built within the `ggplot2` framework

`ggGenomicTracks` provides genomic tracks as standard
[`ggplot2`](https://ggplot2.tidyverse.org) layers that can be stacked,
themed, and combined like any other plot: copy number (total and
allele-specific), structural variant junctions, binned coverage, genes,
mutations, chromatin states, and ideograms. Coordinate systems place one
or more genomic windows, cis or trans, or the whole genome on a single
x-axis, so junctions can be drawn between loci in the style of
[gTrack](https://github.com/mskilab-org/gTrack). It is the graphical
companion to [`devgru`](https://pblaney.github.io/devgru/), and by
default is anchored to hg38.

Full documentation and the function reference are available at
<https://pblaney.github.io/ggGenomicTracks/>.

### Installation

You can install the current development version of `ggGenomicTracks`
with:

``` r

devtools::install_github("pblaney/ggGenomicTracks")
```

Some dependencies are only available on GitHub
([`pblaney/devgru`](https://github.com/pblaney/devgru) and
[`mskilab-org/gUtils`](https://github.com/mskilab-org/gUtils)); these
are declared as `Remotes` and resolved automatically by
`devtools`/`remotes`.

### Function Overview

| Category | Functions |
|----|----|
| **Coordinate systems** | [`coord_locus()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_locus.md) for one window, [`coord_loci()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_loci.md) for several windows on one axis, [`coord_genome()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_genome.md) for the whole genome |
| **Stacking tracks** | [`stack_tracks()`](https://pblaney.github.io/ggGenomicTracks/reference/stack_tracks.md) aligns track plots on shared coordinates; [`highlight_stack()`](https://pblaney.github.io/ggGenomicTracks/reference/highlight_stack.md) marks regions across all tracks |
| **Copy number** | [`geom_copy_number()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_copy_number.md) for total or allele-specific integer copy number |
| **Structural variants** | [`geom_junction()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_junction.md) for cis and trans junction arcs, [`geom_loose_end()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_loose_end.md) for loose ends |
| **Coverage and signal** | [`geom_coverage()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_coverage.md), [`geom_ridge_track()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_ridge_track.md), [`geom_bed_track()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_bed_track.md), [`geom_epilogo_col()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_epilogo_col.md) |
| **Mutations** | [`geom_mutation_count()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_mutation_count.md), [`geom_vaf_landscape()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_vaf_landscape.md) |
| **Genes and annotations** | [`geom_gene_track()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_gene_track.md), [`annotation_locus_scale()`](https://pblaney.github.io/ggGenomicTracks/reference/annotation_locus_scale.md) |
| **Ideograms** | [`geom_ideogram()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_ideogram.md), [`plot_ideogram()`](https://pblaney.github.io/ggGenomicTracks/reference/plot_ideogram.md), [`annotate_ideogram_window()`](https://pblaney.github.io/ggGenomicTracks/reference/annotate_ideogram_window.md), [`annotate_ideogram_zoom()`](https://pblaney.github.io/ggGenomicTracks/reference/annotate_ideogram_zoom.md) |
| **Scales and themes** | [`scale_fill_copynumber()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_colour_copynumber.md), [`scale_colour_sv_class()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_colour_sv_class.md), [`scale_colour_jabba_event()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_colour_jabba_event.md), [`scale_fill_chromhmm()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_fill_chromhmm.md), [`scale_fill_giemsa()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_fill_giemsa.md), [`theme_genomic()`](https://pblaney.github.io/ggGenomicTracks/reference/theme_genomic.md) |
| **Data preparation** | [`junction_copy_number()`](https://pblaney.github.io/ggGenomicTracks/reference/junction_copy_number.md), [`loose_end_copy_number()`](https://pblaney.github.io/ggGenomicTracks/reference/loose_end_copy_number.md), [`annotate_junctions()`](https://pblaney.github.io/ggGenomicTracks/reference/annotate_junctions.md), [`epilogos_long()`](https://pblaney.github.io/ggGenomicTracks/reference/epilogos_long.md) |

See the [reference
page](https://pblaney.github.io/ggGenomicTracks/reference/index.html)
for full details on each function.

### Quick Start

Draw copy number with junctions anchored at the copy number of each
breakpoint, across several windows:

``` r

library(ggGenomicTracks)
library(ggplot2)

cn <- data.table::fread("sample.jabba.cnv.bed")
data.table::setnames(cn, "chrom", "seqnames")
junctions <- devgru::read_bedpe_file("sample.jabba.junctions.bedpe", output_type = "dt")
junctions <- annotate_junctions(junction_copy_number(junctions, cn))

windows <- c(MYC = "chr8:126000001-130000000", IGH = "chr14:105000001-107000000")

ggplot() +
  geom_copy_number(data = cn) +
  geom_junction(data = junctions, aes(y = cn1, yend = cn2, colour = sv_class)) +
  scale_colour_sv_class() +
  coord_loci(windows) +
  theme_genomic()
```

Stack copy number, allele-specific copy number, and coverage tracks on
shared coordinates, here genome-wide:

``` r

coverage <- readRDS("sample.fragcounter.cov.rds")

stack_tracks(
  ggplot() + geom_copy_number(data = cn) + theme_genomic(),
  ggplot() + geom_copy_number(data = cn, alleles = TRUE) + theme_genomic(),
  ggplot(coverage) + geom_coverage(y_range = c(0, 3)) + theme_genomic(),
  windows = "genome",
  heights = c(1.3, 1, 0.8)
)
```

### Bundled Data

`ggGenomicTracks` ships with hg38 cytobands (`cytobands_hg38`) for
ideograms. The raw file it was built from is in `inst/extdata/`.

### Citation

To cite `ggGenomicTracks` in publications, use:

``` r

citation("ggGenomicTracks")
```

### Contributing

Bug reports and feature requests are welcome on the [issue
tracker](https://github.com/pblaney/ggGenomicTracks/issues).
