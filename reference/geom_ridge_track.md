# Signal track

Draws a quantitative signal over genomic ranges, such as ChIP-seq
coverage, a bedGraph, or GISTIC scores, as a shaded area under a step
curve. Each range is drawn flat at its value from `start` to `end`.
Ranges further apart than `max_gap` are not joined, so missing data is
left blank.

There are two modes, set by the mapping:

- **Overlay**: map the signal to `y` to draw curves on a shared
  continuous axis, e.g. `aes(y = score, fill = group)`. The area is
  shaded between 0 and the signal.

- **Ridges**: map a discrete variable to `y` and the signal to `height`
  to draw one curve per row, e.g. `aes(y = sample, height = signal)`.
  Heights are scaled so the largest absolute value spans `ridge_scale`
  rows.

`xmin` and `xmax` default to the `start` and `end` columns. For single
position data, map `x` instead, e.g. `aes(x = position, y = signal)`,
and the points are joined directly.

## Usage

``` r
geom_ridge_track(
  mapping = NULL,
  data = NULL,
  ridge_scale = 0.9,
  max_gap = 1,
  outline.type = "upper",
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

  GRanges object or data frame of ranges with `seqnames`, `start`, and
  `end` columns plus one or more signal columns, in long format with one
  row per range and group.

- ridge_scale:

  In ridge mode, the height of the largest absolute signal in rows,
  default: 0.9. Values above 1 let ridges overlap.

- max_gap:

  Largest distance in bp between the end of one range and the start of
  the next for the two to be joined, default: 1. Use `Inf` to always
  join.

- outline.type:

  Which sides of the area to outline when `colour` is set, one of
  `"upper"`, `"lower"`, `"both"`, or `"full"`. See
  [`ggplot2::geom_ribbon()`](https://ggplot2.tidyverse.org/reference/geom_ribbon.html).

- ...:

  Other arguments passed to
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html),
  e.g. `colour`, `fill`, `alpha`, or `linewidth`.

- na.rm, show.legend, inherit.aes:

  See
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html).

## Value

ggplot2 layer

## Examples

``` r
if (FALSE) { # \dontrun{
# GISTIC scores, deletions drawn below zero; segments are ~10 kb apart so join them
gistic <- data.table::fread("scores.gistic", check.names = TRUE)
gistic <- gistic[, .(seqnames = Chromosome, start = Start, end = End, Type,
                     G.score = ifelse(Type == "Del", -G.score, G.score))]

ggplot2::ggplot(gistic) +
  geom_ridge_track(ggplot2::aes(y = G.score, fill = Type, group = Type), max_gap = 1e5) +
  coord_locus("chr13:48303751-48481890", flank = 1e6) +
  theme_genomic()

# ChIP-seq signal, one ridge per sample
chip <- devgru::read_bed_file("signal.bedGraph", additional_col_names = c("sample", "signal"))

ggplot2::ggplot(chip) +
  geom_ridge_track(ggplot2::aes(y = sample, height = signal, fill = sample)) +
  coord_locus("chr1:203200001-203893000") +
  theme_genomic()
} # }
```
