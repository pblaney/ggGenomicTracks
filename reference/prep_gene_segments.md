# Prepare gene segment data for plotting

Converts to a data frame, subsets to genes of interest, and adds the
strand direction column used by the exon and intron layers. For internal
use only.

## Usage

``` r
prep_gene_segments(data, genes_of_interest = NULL)
```

## Arguments

- data:

  GRanges object or data frame of gene segments with `seqnames`,
  `start`, `end`, `strand`, `gene_name`, and `segment` (`exon`,
  `intron`, `five_prime_utr`, `three_prime_utr`) columns. Defaults to
  [`devgru::gene_body_hg38`](https://pblaney.github.io/devgru/reference/gene_body_hg38.html).

- genes_of_interest:

  Character vector of gene names to draw. Default `NULL` draws every
  gene in the window.
