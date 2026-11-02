# Package index

## Coordinate systems

Place one or more genomic windows, or the whole genome, on a single
x-axis so tracks and junctions can be drawn across loci.

- [`coord_locus()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_locus.md)
  : Cartesian coordinates zoomed to a genomic locus
- [`coord_loci()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_loci.md)
  : Cartesian coordinates across multiple genomic loci
- [`coord_genome()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_genome.md)
  : Cartesian coordinates across the whole genome

## Stacking tracks

Combine track plots that share genomic coordinates into one figure, and
highlight regions across all tracks.

- [`stack_tracks()`](https://pblaney.github.io/ggGenomicTracks/reference/stack_tracks.md)
  : Stack tracks over shared genomic windows
- [`highlight_stack()`](https://pblaney.github.io/ggGenomicTracks/reference/highlight_stack.md)
  : Highlight sub-regions across a stack of tracks

## Layers

### Copy number

- [`geom_copy_number()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_copy_number.md)
  : Copy number segment track

### Structural variants

- [`geom_junction()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_junction.md)
  : Structural variant junction track
- [`geom_loose_end()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_loose_end.md)
  : Loose end track

### Coverage and signal

- [`geom_coverage()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_coverage.md)
  : Binned coverage track
- [`geom_ridge_track()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_ridge_track.md)
  : Signal track
- [`geom_bed_track()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_bed_track.md)
  : Genomic interval track
- [`geom_epilogo_col()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_epilogo_col.md)
  : Epilogos chromatin state track

### Mutations

- [`geom_mutation_count()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_mutation_count.md)
  : Mutation count lollipop track
- [`geom_vaf_landscape()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_vaf_landscape.md)
  : Variant allele frequency track

### Genes and annotations

- [`geom_gene_track()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_gene_track.md)
  : Gene structure track
- [`annotation_locus_scale()`](https://pblaney.github.io/ggGenomicTracks/reference/annotation_locus_scale.md)
  : Locus coordinate and scale bar annotation

## Ideograms

Chromosome ideograms with cytobands, marking the windows shown in a
plot.

- [`geom_ideogram()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_ideogram.md)
  : Chromosome ideogram
- [`plot_ideogram()`](https://pblaney.github.io/ggGenomicTracks/reference/plot_ideogram.md)
  : Ideogram of a locus' chromosome with zoom trace
- [`annotate_ideogram_window()`](https://pblaney.github.io/ggGenomicTracks/reference/annotate_ideogram_window.md)
  : Highlight a window on an ideogram
- [`annotate_ideogram_zoom()`](https://pblaney.github.io/ggGenomicTracks/reference/annotate_ideogram_zoom.md)
  : Zoom trace from an ideogram window to the track below
- [`theme_ideogram()`](https://pblaney.github.io/ggGenomicTracks/reference/theme_ideogram.md)
  : Theme for ideogram plots

## Scales

Color and fill scales for copy number, junction classes, chromatin
states, and cytobands.

- [`scale_colour_copynumber()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_colour_copynumber.md)
  [`scale_fill_copynumber()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_colour_copynumber.md)
  : Integer copy number color scales
- [`scale_colour_cn_state()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_colour_cn_state.md)
  [`scale_fill_cn_state()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_colour_cn_state.md)
  : Copy number state color scales
- [`scale_colour_sv_class()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_colour_sv_class.md)
  : Junction orientation color scale
- [`scale_colour_jabba_event()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_colour_jabba_event.md)
  : JaBbA complex event color scale
- [`scale_colour_sv_type()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_colour_sv_type.md)
  : Structural variant type color scale
- [`scale_fill_chromhmm()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_fill_chromhmm.md)
  : Chromatin state fill scale
- [`scale_fill_giemsa()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_fill_giemsa.md)
  : Giemsa stain fill scale for ideograms

## Themes

- [`theme_genomic()`](https://pblaney.github.io/ggGenomicTracks/reference/theme_genomic.md)
  : Theme for genomic track plots
- [`label_drop_first()`](https://pblaney.github.io/ggGenomicTracks/reference/label_drop_first.md)
  : Blank the lowest axis label

## Data preparation

Helpers that add columns used by the layers, e.g. copy number at each
junction breakpoint.

- [`junction_copy_number()`](https://pblaney.github.io/ggGenomicTracks/reference/junction_copy_number.md)
  : Copy number at each junction breakpoint
- [`loose_end_copy_number()`](https://pblaney.github.io/ggGenomicTracks/reference/loose_end_copy_number.md)
  : Copy number at each loose end
- [`annotate_junctions()`](https://pblaney.github.io/ggGenomicTracks/reference/annotate_junctions.md)
  : Classify junctions by orientation and complex event
- [`epilogos_long()`](https://pblaney.github.io/ggGenomicTracks/reference/epilogos_long.md)
  : Convert epilogos scores to long format
- [`get_linear_ranges_per_chromosome()`](https://pblaney.github.io/ggGenomicTracks/reference/get_linear_ranges_per_chromosome.md)
  : Linearization of GRanges object for single chromosome
- [`gr_absolute_ranges()`](https://pblaney.github.io/ggGenomicTracks/reference/gr_absolute_ranges.md)
  : Linearization of GRanges object

## Data

- [`cytobands_hg38`](https://pblaney.github.io/ggGenomicTracks/reference/cytobands_hg38.md)
  : Cytoband annotations per chromosome on hg38
