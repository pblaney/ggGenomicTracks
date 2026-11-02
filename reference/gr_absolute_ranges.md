# Linearization of GRanges object

Run linearization of GRanges object where the input location relative to
the individual chromosome that is transformed to an absolute location on
a linear end-to-end genome. Adds 2 columns to input GRanges object:

- start_abs_linear

- end_abs_linear

## Usage

``` r
gr_absolute_ranges(
  input_gr,
  cpus = 1,
  seq_lengths = gUtils::hg_seqlengths()[1:24]
)
```

## Arguments

- input_gr:

  GenomicRanges object to transform

- cpus:

  Number of CPUs to use for parallel execution, default: 1

- seq_lengths:

  Named vector object used as template for new seq details.

## Value

Input GRanges object with absolute start and end in linear genome added
