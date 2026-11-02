# Copy number of the segment containing each position

For internal use only.

## Usage

``` r
lookup_copy_number(
  seqnames,
  position,
  segments,
  cn_column = "total_cn",
  cn_cap = NULL
)
```

## Arguments

- seqnames, position:

  Chromosomes and positions to look up

- segments:

  GRanges or data frame of copy number segments with `seqnames`,
  `start`, and `end` columns, for a single sample.

- cn_column:

  Column of `segments` with the copy number, default: `"total_cn"`

- cn_cap:

  Highest copy number, matching `cn_cap` of
  [`geom_copy_number()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_copy_number.md)
  so arcs meet capped segments. Default: no cap.

## Value

Copy number per position, `NA` where no segment contains it
