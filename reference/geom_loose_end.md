# Loose end track

Draws loose ends, breakpoints with copy number change but no junction
partner as inferred by JaBbA, as a short stem topped with an open
circle. `x` defaults to the `start` column. Map `y` to the copy number
at each loose end, e.g. from
[`loose_end_copy_number()`](https://pblaney.github.io/ggGenomicTracks/reference/loose_end_copy_number.md),
to place them on a copy number track; otherwise they sit at the bottom
of the panel.

## Usage

``` r
geom_loose_end(
  mapping = NULL,
  data = NULL,
  stem_height = 3,
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

  Data frame of loose ends with `seqnames` and `start` columns, e.g. a
  JaBbA `loose_ends.bed`.

- stem_height:

  Height of the stems in mm, default: 3

- ...:

  Other arguments passed to
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html),
  e.g. `colour` or `linewidth`.

- na.rm, show.legend, inherit.aes:

  See
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html).

## Value

ggplot2 layer

## Examples

``` r
if (FALSE) { # \dontrun{
loose <- data.table::fread("sample.jabba.phased.loose_ends.bed",
                           col.names = c("seqnames", "start", "end"))

ggplot2::ggplot() +
  geom_copy_number(data = cn) +
  geom_loose_end(data = loose_end_copy_number(loose, cn), ggplot2::aes(y = cn)) +
  coord_genome() +
  theme_genomic()
} # }
```
