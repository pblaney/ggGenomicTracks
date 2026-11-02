# Theme for genomic track plots

Minimal theme for a single genomic track: blank panel, y-axis line with
a horizontal title, and tight margins for stacking tracks. The x-axis
style is set with `x_axis`; use `"full"` on the bottom track of a stack
and `"ticks"` or `"blank"` on the tracks above it. The top axis is left
for the window labels of
[`coord_loci()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_loci.md).
Use `y_axis = "blank"` for tracks labeled within the panel, such as
[`geom_gene_track()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_gene_track.md).

## Usage

``` r
theme_genomic(
  x_axis = c("full", "ticks", "blank"),
  y_axis = c("full", "blank"),
  base_size = 11
)
```

## Arguments

- x_axis:

  Style of the x-axis, one of `"full"`, `"ticks"`, or `"blank"`.

- y_axis:

  Style of the y-axis, `"full"` or `"blank"` (no line, title, text, or
  ticks).

- base_size:

  Base font size, default: 11

## Value

ggplot2 theme object
