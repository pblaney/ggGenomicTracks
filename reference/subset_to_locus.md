# Subset a layer's data to a genomic locus

Keeps rows on the locus chromosome that overlap the locus. Chromosome
names are compared without the `chr` prefix so `1` and `chr1` style data
both match. Data without a `seqnames` column is returned unchanged. For
internal use only.

## Usage

``` r
subset_to_locus(data, locus)
```

## Arguments

- data:

  Data frame of a single layer

- locus:

  List with `chrom`, `start`, and `end` of the locus
