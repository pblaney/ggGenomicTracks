# Cartesian coordinates across multiple genomic loci

Places several genomic windows side by side on one linear x-axis, so
tracks can show distant loci, on the same or different chromosomes, in a
single panel. Every layer whose data has `seqnames` and `start` columns
has its positions converted to the linear axis before aesthetics are
evaluated: rows outside all windows are dropped, and with more than one
window, ranges are clipped to the window edges and split across windows
they span. Data with two anchors, with `seqnames1`, `start1`,
`seqnames2`, and `start2` columns as from
`devgru::read_bedpe_file(output_type = "dt")`, has both anchors
converted and keeps rows with at least one anchor in view; the anchor
out of view is set to `NA`. Layers without these columns are left
untouched.

Also adds an x-axis scale with breaks and labels in each window's own
coordinates, and window labels on the top axis.

## Usage

``` r
coord_loci(
  windows,
  flank = 0,
  widths = c("proportional", "equal"),
  gap = 0.02,
  dividers = TRUE,
  window_labels = NULL,
  expand = TRUE,
  clip = "on"
)
```

## Arguments

- windows:

  Genomic windows as one of:

  - A character vector of windows, e.g. `"chr8:127600001-128000000"`.
    Names, if given, are used as window labels.

  - A data frame with `seqnames`, `start`, and `end` columns and an
    optional `label` column.

  - A GRanges object, with an optional `label` metadata column.

- flank:

  Number of bp to extend both ends of each window, default: 0

- widths:

  How the axis is shared between windows: `"proportional"` (default) to
  each window's length in bp, or `"equal"` so small windows stay
  readable next to large ones, each then at its own scale.

- gap:

  Space between windows as a fraction of the total axis, default: 0.02

- dividers:

  Draw dashed lines between windows, default: TRUE

- window_labels:

  Labels of the windows on the top axis, as one of:

  - `NULL` (default): labeled with more than one window, by the names of
    `windows` or else the chromosome.

  - A character vector with one label per window.

  - `TRUE` or `FALSE` to show or hide the default labels. In a stack,
    keep the labels on the top track only.

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
windows <- c(MYC = "chr8:127700001-128000000", IGH = "chr14:105550001-106880000")

ggplot2::ggplot() +
  geom_gene_track(genes_of_interest = c("MYC", "PVT1")) +
  coord_loci(windows) +
  theme_genomic(y_axis = "blank")
} # }
```
