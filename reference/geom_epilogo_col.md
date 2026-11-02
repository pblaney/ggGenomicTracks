# Epilogos chromatin state track

Draws epilogos per-state scores as stacked bars spanning each bin, one
color per chromatin state. Positive and negative scores, as in epilogos
group comparisons, are stacked separately above and below zero. Data is
in long format with one row per bin and state; use
[`epilogos_long()`](https://pblaney.github.io/ggGenomicTracks/reference/epilogos_long.md)
to convert epilogos output. `xmin`, `xmax`, `y`, and `fill` default to
the `start`, `end`, `score`, and `state` columns.

Pair with
[`scale_fill_chromhmm()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_fill_chromhmm.md)
for state colors. Also adds an "Epilogos" y-axis title; add your own
[`labs()`](https://ggplot2.tidyverse.org/reference/labs.html) after it
to override it.

## Usage

``` r
geom_epilogo_col(
  mapping = NULL,
  data = NULL,
  background = "black",
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

- data:

  GRanges object or data frame of epilogos scores in long format, e.g.
  from
  [`epilogos_long()`](https://pblaney.github.io/ggGenomicTracks/reference/epilogos_long.md).

- background:

  Fill color of the panel behind the bars, default: `"black"`. Use `NA`
  for none.

- ...:

  Other arguments passed to
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html),
  e.g. `alpha`.

- na.rm, show.legend, inherit.aes:

  See
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html).

## Value

List of ggplot2 layers and a y-axis title, added to a plot with `+`

## Examples

``` r
if (FALSE) { # \dontrun{
scores <- data.table::fread("scores_mm_epilogos.chr3.txt.gz")

ggplot2::ggplot(epilogos_long(scores)) +
  geom_epilogo_col() +
  scale_fill_chromhmm() +
  coord_locus("chr3:187734001-187751500", flank = 5000) +
  theme_genomic()
} # }
```
