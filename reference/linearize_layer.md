# Convert a layer's data to the linear axis

Converts `start`/`end` positions, or both anchors of two-anchor data, to
the linear axis of the window layout, and drops rows out of view. With
more than one window, ranges are clipped to window edges and split
across windows. Data without position columns but with `seqnames` is
subset to the window chromosomes. For internal use only.

## Usage

``` r
linearize_layer(data, layout)
```

## Arguments

- data:

  Data frame of a single layer

- layout:

  Window layout from
  [`locus_layout()`](https://pblaney.github.io/ggGenomicTracks/reference/locus_layout.md)
