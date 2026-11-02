# Chromosome ideogram

Draws the cytobands of a single chromosome between y = 0 and y = 1,
colored by Giemsa stain with a pinched centromere and an outline. Pair
with
[`scale_fill_giemsa()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_fill_giemsa.md)
for the standard band colors,
[`annotate_ideogram_window()`](https://pblaney.github.io/ggGenomicTracks/reference/annotate_ideogram_window.md)
and
[`annotate_ideogram_zoom()`](https://pblaney.github.io/ggGenomicTracks/reference/annotate_ideogram_zoom.md)
to mark a locus, and
[`theme_ideogram()`](https://pblaney.github.io/ggGenomicTracks/reference/theme_ideogram.md).
See
[`plot_ideogram()`](https://pblaney.github.io/ggGenomicTracks/reference/plot_ideogram.md)
for all of these in one call.

## Usage

``` r
geom_ideogram(chrom, data = NULL, outline_colour = "gray20", linewidth = 0.3)
```

## Arguments

- chrom:

  Chromosome to draw, e.g. `"chr3"`.

- data:

  GRanges object or data frame of cytobands with `seqnames`, `start`,
  `end`, `name`, and `gieStain` columns. Defaults to
  [cytobands_hg38](https://pblaney.github.io/ggGenomicTracks/reference/cytobands_hg38.md).

- outline_colour:

  Color of the chromosome outline.

- linewidth:

  Width of the chromosome outline.

## Value

List of ggplot2 layers

## Examples

``` r
if (FALSE) { # \dontrun{
ggplot2::ggplot() +
  geom_ideogram("chr3") +
  annotate_ideogram_window("chr3:187734001-187751500") +
  scale_fill_giemsa() +
  theme_ideogram()
} # }
```
