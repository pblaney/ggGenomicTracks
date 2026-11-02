# Copy number at each junction breakpoint

Adds the copy number on the joined side of each breakpoint, `cn1` and
`cn2`, for anchoring
[`geom_junction()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_junction.md)
arcs on a copy number track. A breakpoint on the `+` strand joins the
segment ending at it, on its left, and one on the `-` strand joins the
segment starting after it, on its right; without strand columns the
segment containing the breakpoint is used.

## Usage

``` r
junction_copy_number(
  junctions,
  segments,
  cn_column = "total_cn",
  cn_cap = NULL
)
```

## Arguments

- junctions:

  Data frame of junctions with `seqnames1`, `start1`, `seqnames2`, and
  `start2` columns, and optionally `strand1` and `strand2`.

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

`junctions` as a data frame with `cn1` and `cn2` columns added
