# Highlight a window on an ideogram

Draws a box over the window on a
[`geom_ideogram()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_ideogram.md).
The outline keeps the box visible even when the window is far narrower
than a pixel at chromosome scale.

## Usage

``` r
annotate_ideogram_window(
  window_of_interest,
  flank = 0,
  fill = "#FFA319",
  colour = fill,
  linewidth = 1
)
```

## Arguments

- window_of_interest:

  Genomic window as a string, e.g. `"chr3:187734001-187751500"`.

- flank:

  Number of bp to extend both ends of the window, default: 0

- fill:

  Fill color of the box.

- colour:

  Outline color of the box.

- linewidth:

  Outline width of the box.

## Value

ggplot2 layer
