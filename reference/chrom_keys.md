# Chromosome names without the chr prefix

Strips the prefix from the unique names only, so it stays fast for
millions of rows. For internal use only.

## Usage

``` r
chrom_keys(seqnames)
```

## Arguments

- seqnames:

  Chromosome names, as a character vector or factor
