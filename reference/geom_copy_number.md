# Copy number segment track

Draws copy number segments as thin horizontal bars at their integer copy
number, as in gTrack. `xmin`, `xmax`, and `y` default to the `start`,
`end`, and `total_cn` columns, as in JaBbA copy number profiles; map `y`
to another column such as `major_allele` to change it. Segments are
filled by their integer copy number with
[`scale_fill_copynumber()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_colour_copynumber.md)
unless `fill` is mapped or set, e.g. by state with
`aes(fill = cn_state)` and
[`scale_fill_cn_state()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_colour_cn_state.md).
Bars use `fill` so `colour` stays free for
[`geom_junction()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_junction.md)
in the same plot.

Segments above `cn_cap` are drawn at the cap and marked with a triangle,
so a few highly amplified segments do not flatten the rest of the
profile. Adjacent segments are joined by thin gray reference connectors,
as in JaBbA genome graphs.

With `alleles`, the major and minor allele copy numbers are drawn as two
semi-transparent series, red and blue, so where they are equal the bars
overlap in purple.

Also adds a y-axis scale with whole-number breaks starting at zero and a
"Copy Number" y-axis title. Add your own
[`scale_y_continuous()`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)
or [`labs()`](https://ggplot2.tidyverse.org/reference/labs.html) after
it to override them.

## Usage

``` r
geom_copy_number(
  mapping = NULL,
  data = NULL,
  cn_cap = NULL,
  bar_height = 1.5,
  connectors = TRUE,
  alleles = NULL,
  show_first_y_label = TRUE,
  ...,
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
)
```

## Arguments

- mapping:

  Set of aesthetic mappings created by
  [`ggplot2::aes()`](https://ggplot2.tidyverse.org/reference/aes.html),
  e.g. `fill`.

- data:

  GRanges object or data frame of copy number segments with `seqnames`,
  `start`, `end`, and `total_cn` columns, e.g. a JaBbA `cnv.bed`
  profile.

- cn_cap:

  Highest copy number shown, default: no cap. Segments above it are
  drawn at the cap with a triangle marker.

- bar_height:

  Height of the bars in mm, default: 1.5

- connectors:

  Join adjacent segments with thin gray lines, default: TRUE. Segments
  are assumed to come from one sample per layer.

- alleles:

  Draw allele-specific copy number instead of total copy number, by
  giving the major and minor allele columns, e.g.
  `c("major_allele", "minor_allele")`, or `TRUE` for those JaBbA column
  names. Default: `NULL` for total copy number.

- show_first_y_label:

  Label the lowest y-axis break, default: TRUE. See
  [`label_drop_first()`](https://pblaney.github.io/ggGenomicTracks/reference/label_drop_first.md).

- ...:

  Other arguments passed to
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html),
  e.g. `fill`, `colour` (bar outline), or `alpha`.

- na.rm, show.legend, inherit.aes:

  See
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html).

## Value

List of a ggplot2 layer, a y-axis scale, a y-axis title, and a copy
number fill scale unless `fill` is given, added to a plot with `+`

## Examples

``` r
if (FALSE) { # \dontrun{
cn <- data.table::fread("sample.jabba.phased.cnv.bed")
data.table::setnames(cn, "chrom", "seqnames")

ggplot2::ggplot(cn) +
  geom_copy_number(cn_cap = 8) +
  coord_genome() +
  theme_genomic()
} # }
```
