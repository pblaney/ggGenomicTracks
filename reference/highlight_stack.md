# Highlight sub-regions across a stack of tracks

Shades one or more sub-regions as unbroken bands running through every
track of a stack, including the space between tracks. Bands cover each
run of adjacent tracks that share the same
[`coord_locus()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_locus.md)
or
[`coord_loci()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_loci.md)
windows; other plots in the stack, such as the ideogram from
[`plot_ideogram()`](https://pblaney.github.io/ggGenomicTracks/reference/plot_ideogram.md),
are left alone. Apply it as the last step, after the stack is assembled
with patchwork.

Bands are drawn over the data by default, so they stay visible on tracks
whose data covers the whole panel, like dense bed tracks. Adjust `alpha`
to set how strongly they tint the data, or use `layer = "below"` to keep
the data untinted.

## Usage

``` r
highlight_stack(
  stack,
  regions,
  fill = NULL,
  alpha = 0.3,
  colour = NA,
  layer = c("above", "below"),
  tracks = NULL,
  show_legend = NULL,
  legend_title = "Highlight"
)
```

## Arguments

- stack:

  A patchwork of tracks, e.g. from
  `patchwork::wrap_plots(..., ncol = 1)`, or a single ggplot.

- regions:

  Regions to highlight as one of:

  - A character vector of windows, e.g. `"chr3:187744872-187744917"`.
    Names, if given, are used as labels.

  - A data frame with `start` and `end` columns and optional `seqnames`,
    `label`, and `fill` columns.

  - A GRanges object, with optional `label` and `fill` metadata columns.

- fill:

  Fill colors, either one color for all regions or one per label (named
  or in order of appearance). Ignored if `regions` has a `fill` column.
  By default a single region is orange and labeled regions use a hue
  palette.

- alpha:

  Opacity of the bands, default: 0.3

- colour:

  Outline color of the bands, default: none

- layer:

  Draw the bands `"above"` (default) or `"below"` the track data.

- tracks:

  Positions of the plots in the stack to cover, e.g. `3:5`. By default
  all plots using
  [`coord_locus()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_locus.md)
  or
  [`coord_loci()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_loci.md).
  Bands still run unbroken through adjacent tracks.

- show_legend:

  Add a legend of the region labels, default: TRUE when regions are
  labeled.

- legend_title:

  Title of the legend, default: "Highlight"

## Value

A gtable of the highlighted stack with class `highlighted_stack`. Print
it to draw, save it with
[`ggplot2::ggsave()`](https://ggplot2.tidyverse.org/reference/ggsave.html),
or wrap it with `patchwork::wrap_elements(full = )` to combine with
other plots.

## Examples

``` r
if (FALSE) { # \dontrun{
win <- "chr3:187744850-187745350"

stack <- patchwork::wrap_plots(
  plot_ideogram(win, flank = 250),
  ggplot2::ggplot(snvs) + geom_mutation_count() + coord_locus(win, 250) +
    theme_genomic(x_axis = "ticks"),
  ggplot2::ggplot() + geom_gene_track(genes_of_interest = "BCL6") +
    coord_locus(win, 250) + theme_genomic(y_axis = "blank"),
  ncol = 1, guides = "collect")

highlight_stack(stack,
                c("Normal MBC Hotspot 1" = "chr3:187744872-187744917",
                  "Normal MBC Hotspot 2" = "chr3:187745199-187745246",
                  "DLBCL Hotspot" = "chr3:187744678-187744698"),
                fill = c("#9BD08A", "#C2EFB4", "#DC766F"),
                legend_title = "Windows")
} # }
```
