# Allele-specific copy number in long format

Stacks the major and minor allele columns into one `allele_cn` column
with an `allele` column of `"Major"` or `"Minor"`. For internal use
only.

## Usage

``` r
allele_long(data, alleles)
```

## Arguments

- data:

  Data frame or GRanges of copy number segments

- alleles:

  Names of the major and minor allele columns
