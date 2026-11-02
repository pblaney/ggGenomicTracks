# Mutation count lollipop track

Counts mutations at each position and draws the counts as lollipops: a
stem from zero and a point sized by the count. `x` defaults to the
`start` column. Counts are made within each group, so mapping `fill` or
`group` (e.g. to `Reference_Allele` or a sample cohort) counts each
level separately.

Also adds a y-axis scale with whole-number breaks starting at zero and a
"Mutation Count" y-axis title. Add your own
[`scale_y_continuous()`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)
or [`labs()`](https://ggplot2.tidyverse.org/reference/labs.html) after
it to override them.

## Usage

``` r
geom_mutation_count(
  mapping = NULL,
  data = NULL,
  size_range = c(1.5, 4),
  show_first_y_label = TRUE,
  stem_colour = "gray10",
  stem_alpha = 0.65,
  stem_linewidth = 0.5,
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

  GRanges object or data frame of mutations with `seqnames` and `start`
  columns, one row per mutation, e.g. from a MAF or VCF read with
  devgru.

- size_range:

  Range of point sizes in mm for the smallest and largest counts,
  default: `c(1.5, 4)`. Point area grows with the count. Map `size` to
  `after_stat(count)` instead to get a size legend.

- show_first_y_label:

  Label the lowest y-axis break, default: TRUE. Set to `FALSE` to blank
  it, keeping the tick, so it does not collide with the top label of the
  track below in a tightly stacked plot. See
  [`label_drop_first()`](https://pblaney.github.io/ggGenomicTracks/reference/label_drop_first.md).

- stem_colour, stem_alpha, stem_linewidth:

  Color, opacity, and width in mm of the stems.

- ...:

  Other arguments passed to
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html),
  e.g. `fill`, `colour` (point outline), `shape`, or `alpha`.

- na.rm, show.legend, inherit.aes:

  See
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html).

## Value

List of a ggplot2 layer, a y-axis scale, and a y-axis title, added to a
plot with `+`

## Computed variables

- `count`: number of mutations at the position within the group.

## Examples

``` r
if (FALSE) { # \dontrun{
ggplot2::ggplot(snvs) +
  geom_mutation_count(fill = "dodgerblue") +
  coord_locus("chr3:187734001-187751500", flank = 5000) +
  theme_genomic()

# Count each cohort separately
ggplot2::ggplot(snvs) +
  geom_mutation_count(ggplot2::aes(fill = cohort), alpha = 0.8) +
  coord_locus("chr3:187734001-187751500", flank = 5000) +
  theme_genomic()
} # }
```
