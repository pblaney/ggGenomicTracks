# Binned coverage track

Draws read coverage per genomic bin as points, as in gTrack, e.g.
fragCounter or dryclean output. `x` defaults to the middle of each bin
and `y` to the `reads.corrected` column. Whole-genome coverage has
millions of bins, more than can be seen or drawn quickly, so bins are
summarized into `n_bins` bins across the x-axis when there are more than
that; zoomed in, every bin is drawn.

Also adds a "Coverage" y-axis title, and with `y_range`, a y-axis scale.
Add your own
[`labs()`](https://ggplot2.tidyverse.org/reference/labs.html) or scale
after it to override them.

## Usage

``` r
geom_coverage(
  mapping = NULL,
  data = NULL,
  n_bins = 3000,
  summary = c("mean", "median"),
  y_range = NULL,
  point_size = 0.4,
  ...,
  na.rm = TRUE,
  show.legend = NA,
  inherit.aes = TRUE
)
```

## Arguments

- mapping:

  Set of aesthetic mappings created by
  [`ggplot2::aes()`](https://ggplot2.tidyverse.org/reference/aes.html),
  e.g. `y` for another coverage column or `colour`.

- data:

  GRanges object or data frame of bins with `seqnames`, `start`, `end`,
  and coverage columns.

- n_bins:

  Largest number of points across the x-axis, default: 3000. Use `Inf`
  to draw every bin.

- summary:

  Summary of the coverage of bins combined into one point, `"mean"`
  (default) or `"median"`.

- y_range:

  Limits of the y-axis, with coverage outside drawn at the limits.
  Default: the range of the data.

- point_size:

  Size of the points in mm, default: 0.4

- ...:

  Other arguments passed to
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html),
  e.g. `colour` or `alpha`.

- na.rm, show.legend, inherit.aes:

  See
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html).

## Value

List of a ggplot2 layer, a y-axis title, and with `y_range`, a y-axis
scale, added to a plot with `+`

## Examples

``` r
if (FALSE) { # \dontrun{
coverage <- readRDS("sample.fragcounter.cov.rds")

ggplot2::ggplot(coverage) +
  geom_coverage(y_range = c(0, 3)) +
  coord_genome() +
  theme_genomic()
} # }
```
