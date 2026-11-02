# Zoom trace from an ideogram window to the track below

Draws a funnel from the window on a
[`geom_ideogram()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_ideogram.md)
down to the full chromosome width at the bottom of the panel. When the
ideogram is stacked above locus tracks with patchwork, the funnel
connects the window to the zoomed tracks.

## Usage

``` r
annotate_ideogram_zoom(
  window_of_interest,
  flank = 0,
  data = NULL,
  fill = "#FFA319",
  alpha = 0.35
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

- fill:

  Fill color of the funnel.

- alpha:

  Transparency of the funnel.

## Value

ggplot2 layer
