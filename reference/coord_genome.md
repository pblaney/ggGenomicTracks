# Cartesian coordinates across the whole genome

Places whole chromosomes end to end on one linear x-axis, for
genome-wide tracks such as copy number. Works like
[`coord_loci()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_loci.md)
with every chromosome as a window: the x-axis is labeled with chromosome
names, and alternate chromosomes are shaded.

## Usage

``` r
coord_genome(
  chromosomes = c(1:22, "X", "Y"),
  seqlengths = NULL,
  gap = 0,
  dividers = TRUE,
  shade = "gray95",
  expand = FALSE,
  clip = "on"
)
```

## Arguments

- chromosomes:

  Chromosomes to show, in order, with or without the `chr` prefix.
  Default: 1-22, X, and Y.

- seqlengths:

  Named vector of chromosome lengths, default: hg38 from
  BSgenome.Hsapiens.UCSC.hg38.

- gap:

  Space between chromosomes as a fraction of the total axis, default: 0

- dividers:

  Draw dashed lines between chromosomes, default: TRUE

- shade:

  Fill color of every other chromosome, default: `"gray95"`. Use `NA`
  for none.

- expand:

  If `TRUE`, the default, adds a small expansion to the window limits.

- clip:

  Should drawing be clipped to the extent of the plot panel? Default:
  `"on"`

## Value

List of a coord and an x-axis scale, added to a plot with `+`

## Examples

``` r
if (FALSE) { # \dontrun{
ggplot2::ggplot(cn_segments) +
  geom_bed_track(ggplot2::aes(y = sample, fill = cn)) +
  coord_genome() +
  theme_genomic()
} # }
```
