# Chromosome and cytoband span label for a locus

Builds a label like "Chr 3 (q27.3)" or "Chr 3 (q27.2-q27.3)" from the
cytobands overlapping the window. For internal use only.

## Usage

``` r
locus_cytoband_label(window_of_interest, flank = 0, data = NULL)
```

## Arguments

- window_of_interest:

  Genomic window as a string, e.g. `"chr3:187734001-187751500"`.

- flank:

  Number of bp to extend both ends of the window, default: 0

- data:

  GRanges object or data frame of cytobands with `seqnames`, `start`,
  `end`, `name`, and `gieStain` columns. Defaults to
  [cytobands_hg38](https://pblaney.github.io/ggGenomicTracks/reference/cytobands_hg38.md).
