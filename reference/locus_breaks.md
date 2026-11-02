# Axis breaks within each window

Pretty breaks inside each window, in proportion to its share of the
axis, labeled in Mbp for windows wider than 10 kbp and bp otherwise.
With several windows, breaks near window edges are dropped so labels of
neighboring windows do not collide. For internal use only.

## Usage

``` r
locus_breaks(layout, n = 8)
```

## Arguments

- layout:

  Window layout from
  [`locus_layout()`](https://pblaney.github.io/ggGenomicTracks/reference/locus_layout.md)

- n:

  Approximate number of breaks across the whole axis, default: 8, or 5
  for a single window

## Value

List of break positions on the linear axis (`lin`) and their labels
