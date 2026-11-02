# Copy number at each loose end

Adds a `cn` column with the copy number of the segment containing each
loose end, for placing
[`geom_loose_end()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_loose_end.md)
on a copy number track.

## Usage

``` r
loose_end_copy_number(
  loose_ends,
  segments,
  cn_column = "total_cn",
  cn_cap = NULL
)
```

## Arguments

- loose_ends:

  Data frame of loose ends with `seqnames` and `start` columns.

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

`loose_ends` as a data frame with a `cn` column added
