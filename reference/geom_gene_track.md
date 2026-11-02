# Gene structure track

Draws gene models with exons, UTRs, and introns with strand direction
arrows, plus one italic label per gene at its 5' or 3' end. One row per
level of the `y` aesthetic, which defaults to `gene_name`. Pair with
[`coord_locus()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_locus.md)
to zoom to a window and `theme_genomic(y_axis = "blank")` for a clean
gene track.

Labels for genes extending past the window are kept at the panel edge.

## Usage

``` r
geom_gene_track(
  mapping = NULL,
  data = NULL,
  genes_of_interest = NULL,
  exon_fill_color = "darkred",
  gene_name_side = "five_prime",
  gene_label_size = 4,
  gene_label_padding = 1.5,
  exon_height = 7,
  exon_arrowhead_width = 7,
  utr_linewidth = 5,
  intron_linewidth = 0.7,
  inherit.aes = FALSE
)
```

## Arguments

- mapping:

  Set of aesthetic mappings created by
  [`ggplot2::aes()`](https://ggplot2.tidyverse.org/reference/aes.html).
  Defaults to `y = gene_name`. Map `fill` (e.g. `fill = gene_name`) to
  color exons by a column.

- data:

  GRanges object or data frame of gene segments with `seqnames`,
  `start`, `end`, `strand`, `gene_name`, and `segment` (`exon`,
  `intron`, `five_prime_utr`, `three_prime_utr`) columns. Defaults to
  [`devgru::gene_body_hg38`](https://pblaney.github.io/devgru/reference/gene_body_hg38.html).

- genes_of_interest:

  Character vector of gene names to draw. Default `NULL` draws every
  gene in the window.

- exon_fill_color:

  Exon fill color, used when `fill` is not mapped.

- gene_name_side:

  Side to place gene labels, `"five_prime"` or `"three_prime"`. Either a
  single value or one per gene in `genes_of_interest`.

- gene_label_size:

  Gene label text size, default: 4

- gene_label_padding:

  Space between gene end and label in mm, default: 1.5

- exon_height:

  Height of exons in mm, default: 7

- exon_arrowhead_width:

  Width of exon arrowheads in mm. Set to 0 for plain rectangular exons,
  default: 7

- utr_linewidth:

  Line width of UTRs, default: 5

- intron_linewidth:

  Line width of introns, default: 0.7

- inherit.aes:

  See
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html).

## Value

List of ggplot2 layers

## Examples

``` r
if (FALSE) { # \dontrun{
ggplot2::ggplot() +
  geom_gene_track(genes_of_interest = "BCL6") +
  coord_locus("chr3:187734001-187751500", flank = 5000) +
  theme_genomic(y_axis = "blank")

# All genes in the window, exons colored per gene
ggplot2::ggplot() +
  geom_gene_track(ggplot2::aes(fill = gene_name)) +
  coord_locus("chr3:187600001-187800000") +
  theme_genomic(y_axis = "blank")
} # }
```
