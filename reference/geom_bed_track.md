# Genomic interval track

Draws BED-like intervals as bars of fixed height, one row per level of
the `y` aesthetic. `xmin` and `xmax` default to the `start` and `end`
columns. Pair with
[`coord_locus()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_locus.md)
to zoom to a window and
[`theme_genomic()`](https://pblaney.github.io/ggGenomicTracks/reference/theme_genomic.md)
for track styling.

## Usage

``` r
geom_bed_track(
  mapping = NULL,
  data = NULL,
  bar_height = 5,
  ...,
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
)
```

## Arguments

- mapping:

  Set of aesthetic mappings created by
  [`ggplot2::aes()`](https://ggplot2.tidyverse.org/reference/aes.html).
  Typically `y` (e.g. `sample`) and optionally `fill`.

- data:

  GRanges object or data frame of intervals with `seqnames`, `start`,
  and `end` columns, e.g. from
  [`devgru::read_bed_file()`](https://pblaney.github.io/devgru/reference/read_bed_file.html).

- bar_height:

  Height of the bars in mm, default: 5

- ...:

  Other arguments passed to
  [`gggenes::geom_gene_arrow()`](https://wilkox.org/gggenes/reference/geom_gene_arrow.html),
  e.g. `colour` or `alpha`. Bars have no outline unless `colour` is
  mapped or set.

- na.rm, show.legend, inherit.aes:

  See
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html).

## Value

ggplot2 layer

## Examples

``` r
if (FALSE) { # \dontrun{
chrom_states <- devgru::read_bed_file("segments.bed", additional_col_names = "segment_label")

ggplot2::ggplot() +
  geom_bed_track(data = chrom_states, ggplot2::aes(y = sample, fill = segment_label)) +
  coord_locus("chr1:203296001-203313000") +
  theme_genomic()

# Tornado sort: order rows by interval width, largest at top
ggplot2::ggplot() +
  geom_bed_track(data = chrom_states, ggplot2::aes(y = reorder(sample, end - start))) +
  coord_locus("chr1:203296001-203313000")
} # }
```
