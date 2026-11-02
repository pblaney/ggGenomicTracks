# Convert genomic positions in a window to the linear axis

For internal use only.

## Usage

``` r
to_linear(position, layout, window)
```

## Arguments

- position:

  Genomic positions

- layout:

  Window layout from
  [`locus_layout()`](https://pblaney.github.io/ggGenomicTracks/reference/locus_layout.md)

- window:

  Row of the layout the positions belong to, one or one per position
