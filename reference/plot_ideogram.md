# Ideogram of a locus' chromosome with zoom trace

Complete ideogram plot for the top of a track stack: the chromosome with
the window highlighted, a zoom trace down to the tracks below,
chromosome positions on top, and a title of the chromosome and cytoband
span, e.g. "Chr 3 (q27.3)". Built from
[`geom_ideogram()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_ideogram.md)
and its annotation layers, so it can be modified with `+`.

## Usage

``` r
plot_ideogram(
  window_of_interest,
  flank = 0,
  data = NULL,
  highlight_fill = "#FFA319",
  zoom_trace = TRUE
)
```

## Arguments

- window_of_interest:

  Genomic window as a string, e.g. `"chr3:187734001-187751500"`.

- flank:

  Number of bp to extend both ends of the window, default: 0

- data:

  GRanges object or data frame of cytobands with `seqnames`, `start`,
  `end`, `name`, and `gieStain` columns. Defaults to
  [cytobands_hg38](https://pblaney.github.io/ggGenomicTracks/reference/cytobands_hg38.md).

- highlight_fill:

  Color of the window highlight and zoom trace.

- zoom_trace:

  Draw the zoom trace, default: TRUE

## Value

ggplot object

## Examples

``` r
if (FALSE) { # \dontrun{
win <- "chr3:187734001-187751500"
patchwork::wrap_plots(plot_ideogram(win, flank = 5000),
                      bed_plot + coord_locus(win, flank = 5000),
                      ncol = 1, heights = c(0.4, 1))
} # }
```
