# Standardize genomic windows

Converts window strings, data frames, or GRanges to a data frame with
`seqnames`, `start`, `end`, and `label` columns, extended by the flank.
For internal use only.

## Usage

``` r
as_locus_windows(windows, flank = 0)
```

## Arguments

- windows:

  Genomic windows as one of:

  - A character vector of windows, e.g. `"chr8:127600001-128000000"`.
    Names, if given, are used as window labels.

  - A data frame with `seqnames`, `start`, and `end` columns and an
    optional `label` column.

  - A GRanges object, with an optional `label` metadata column.

- flank:

  Number of bp to extend both ends of each window, default: 0
