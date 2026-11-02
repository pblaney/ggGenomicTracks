# Parse a genomic window string

Parses a window string and extends it by the flank on both ends. For
internal use only.

## Usage

``` r
parse_locus(window_of_interest, flank = 0)
```

## Arguments

- window_of_interest:

  Genomic window as a string, e.g. `"chr3:187734001-187751500"`.

- flank:

  Number of bp to extend both ends of the window, default: 0

## Value

List with `chrom`, `start`, and `end` of the locus
