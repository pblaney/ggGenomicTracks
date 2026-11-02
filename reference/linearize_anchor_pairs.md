# Convert both anchors of two-anchor data to the linear axis

Each anchor is converted using the first window it overlaps, or set to
`NA` when it is out of view. Rows with neither anchor in view are
dropped. For internal use only.

## Usage

``` r
linearize_anchor_pairs(data, layout)
```

## Arguments

- data:

  Data frame of a single layer

- layout:

  Window layout from
  [`locus_layout()`](https://pblaney.github.io/ggGenomicTracks/reference/locus_layout.md)
