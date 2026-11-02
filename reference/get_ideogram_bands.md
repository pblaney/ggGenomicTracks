# Cytobands of a single chromosome

Subsets cytoband data to one chromosome, sorted by position. For
internal use only.

## Usage

``` r
get_ideogram_bands(chrom, data = NULL)
```

## Arguments

- chrom:

  Chromosome to draw, e.g. `"chr3"`.

- data:

  GRanges object or data frame of cytobands with `seqnames`, `start`,
  `end`, `name`, and `gieStain` columns. Defaults to
  [cytobands_hg38](https://pblaney.github.io/ggGenomicTracks/reference/cytobands_hg38.md).
