# Cartesian coordinates zoomed to a genomic locus

Sets the x-axis of a plot to a genomic window. Every layer whose data
has a `seqnames` column is subset to rows on the window's chromosome
that overlap the window, before aesthetics are evaluated, so discrete
scales only show levels present in the window. Layers without a
`seqnames` column are left untouched.

Also adds an x-axis scale named after the chromosome (e.g. "Chr 3
(Mbp)") with labels in Mbp for windows wider than 10 kbp and bp
otherwise. Add your own
[`scale_x_continuous()`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)
after `coord_locus()` to override it.

`coord_locus()` is
[`coord_loci()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_loci.md)
with a single window, and x positions are the genomic coordinates
themselves.

## Usage

``` r
coord_locus(window_of_interest, flank = 0, expand = TRUE, clip = "on")
```

## Arguments

- window_of_interest:

  Genomic window as a string, e.g. `"chr3:187734001-187751500"`.

- flank:

  Number of bp to extend both ends of the window, default: 0

- expand:

  If `TRUE`, the default, adds a small expansion to the window limits.

- clip:

  Should drawing be clipped to the extent of the plot panel? Default:
  `"on"`

## Value

List of a coord and an x-axis scale, added to a plot with `+`
