# Convert genomic data to a plain data frame for a ggplot2 layer

Converts GRanges objects to a data frame with `seqnames`, `start`,
`end`, `strand`, and `width` columns plus metadata columns. Data frames
and `NULL` are returned unchanged. For internal use only.

## Usage

``` r
as_genomic_df(data)
```

## Arguments

- data:

  GRanges object, data frame, or `NULL`
