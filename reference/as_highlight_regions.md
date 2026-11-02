# Standardize regions to highlight

Converts window strings, data frames, or GRanges to a data frame with
`seqnames`, `start`, `end`, `label`, and `fill` columns. For internal
use only.

## Usage

``` r
as_highlight_regions(regions, fill = NULL)
```

## Arguments

- regions:

  Regions to highlight as one of:

  - A character vector of windows, e.g. `"chr3:187744872-187744917"`.
    Names, if given, are used as labels.

  - A data frame with `start` and `end` columns and optional `seqnames`,
    `label`, and `fill` columns.

  - A GRanges object, with optional `label` and `fill` metadata columns.

- fill:

  Fill colors, either one color for all regions or one per label (named
  or in order of appearance). Ignored if `regions` has a `fill` column.
  By default a single region is orange and labeled regions use a hue
  palette.
