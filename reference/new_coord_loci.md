# Build the multi-locus coord

Extends
[`ggplot2::coord_cartesian()`](https://ggplot2.tidyverse.org/reference/coord_cartesian.html)
with a window layout used to convert layer data to the linear axis, and
draws window dividers and shading. For internal use only.

## Usage

``` r
new_coord_loci(
  layout,
  expand = TRUE,
  clip = "on",
  dividers = FALSE,
  shade = NA,
  spec = NULL
)
```

## Arguments

- layout:

  Window layout from
  [`locus_layout()`](https://pblaney.github.io/ggGenomicTracks/reference/locus_layout.md)

- expand:

  If `TRUE`, the default, adds a small expansion to the window limits.

- clip:

  Should drawing be clipped to the extent of the plot panel? Default:
  `"on"`

- dividers:

  Draw dashed lines between windows

- shade:

  Fill color of every other window, or `NA` for none

- spec:

  List of the coord function name and its arguments, used by
  [`stack_tracks()`](https://pblaney.github.io/ggGenomicTracks/reference/stack_tracks.md)
  to rebuild the coord for each track
