# Variant allele frequency track

Draws each mutation as a point at its position and variant allele
frequency (VAF), colored by the local density of points so clusters of
mutations at similar VAF stand out. `x` and `y` default to the `start`
and `VAF` columns. Density is estimated from the mutations within the
[`coord_locus()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_locus.md)
window.

Also adds a "Density" color scale labeled from low to high, a y-axis
fixed to 0-1, and a "VAF" y-axis title. Add your own scales or
[`labs()`](https://ggplot2.tidyverse.org/reference/labs.html) after it
to override them.

## Usage

``` r
geom_vaf_landscape(
  mapping = NULL,
  data = NULL,
  point_size = 3,
  background = "gray50",
  method = "kde2d",
  adjust = 1,
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
  [`ggplot2::aes()`](https://ggplot2.tidyverse.org/reference/aes.html).

- data:

  GRanges object or data frame of mutations with `seqnames`, `start`,
  and `VAF` columns, e.g. a MAF with `VAF = t_alt_count / t_depth`.

- point_size:

  Size of the points in mm, default: 3

- background:

  Fill color of the panel behind the points, default: `"gray50"`, which
  gives the light end of the density colors contrast. Use `NA` for none.

- method, adjust:

  Density estimation method and bandwidth adjustment, see
  [`ggpointdensity::stat_pointdensity()`](https://rdrr.io/pkg/ggpointdensity/man/stat_pointdensity.html).
  Default: `"kde2d"`, 1

- show_first_y_label:

  Label the lowest y-axis break, default: TRUE. See
  [`label_drop_first()`](https://pblaney.github.io/ggGenomicTracks/reference/label_drop_first.md).

- ...:

  Other arguments passed to
  [`ggpointdensity::stat_pointdensity()`](https://rdrr.io/pkg/ggpointdensity/man/stat_pointdensity.html),
  e.g. `alpha` or `shape`.

- na.rm, show.legend, inherit.aes:

  See
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html).

## Value

List of ggplot2 layers, scales, and a y-axis title, added to a plot with
`+`

## Examples

``` r
if (FALSE) { # \dontrun{
snvs$VAF <- snvs$t_alt_count / snvs$t_depth

ggplot2::ggplot(snvs) +
  geom_vaf_landscape() +
  coord_locus("chr1:203296001-203313000") +
  theme_genomic()
} # }
```
