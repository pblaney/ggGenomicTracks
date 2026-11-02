# Linearization of GRanges object for single chromosome

Run linearization of GRanges object where the input location relative to
the individual chromosome that is transformed to an absolute location on
a linear end-to-end genome. Useful for graphing across the whole genome.

Note, this function was designed to be run as part of the
[`gr_absolute_ranges()`](https://pblaney.github.io/ggGenomicTracks/reference/gr_absolute_ranges.md)
workflow. Also, see `chrompar()` for executing in parallel.

## Usage

``` r
get_linear_ranges_per_chromosome(
  chrom_for_linearizing,
  input_ranges,
  seq_lengths
)
```

## Arguments

- chrom_for_linearizing:

  The chromosome to run linearization

- input_ranges:

  GenomicRanges object of interest

- seq_lengths:

  Named vector object used as template for new seq details.

## Value

data.table with absolute start and end in linear genome
