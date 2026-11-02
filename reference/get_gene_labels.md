# One label position per gene

Collapses gene segments to one row per gene spanning the whole gene,
with the label anchored at the 5' or 3' end and justified away from the
gene. For internal use only.

## Usage

``` r
get_gene_labels(
  gene_dt,
  genes_of_interest = NULL,
  gene_name_side = "five_prime"
)
```

## Arguments

- gene_dt:

  Data frame from
  [`prep_gene_segments()`](https://pblaney.github.io/ggGenomicTracks/reference/prep_gene_segments.md)

- genes_of_interest:

  Character vector of gene names to draw. Default `NULL` draws every
  gene in the window.

- gene_name_side:

  Side to place gene labels, `"five_prime"` or `"three_prime"`. Either a
  single value or one per gene in `genes_of_interest`.
