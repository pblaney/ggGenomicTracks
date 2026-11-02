# Stack tracks over shared genomic windows

Stacks tracks vertically with patchwork, giving every track the same
genomic windows so their x-axes line up. Each track gets the same
[`coord_loci()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_loci.md)
or
[`coord_genome()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_genome.md),
window labels are kept on the top track only, and the x-axis labels on
the bottom track only. Tracks with another coord, such as the ideogram
from
[`plot_ideogram()`](https://pblaney.github.io/ggGenomicTracks/reference/plot_ideogram.md),
are stacked as they are.

## Usage

``` r
stack_tracks(
  ...,
  windows = NULL,
  flank = 0,
  widths = c("proportional", "equal"),
  gap = 0.02,
  dividers = TRUE,
  window_labels = NULL,
  heights = NULL,
  upper_x_axis = c("ticks", "blank", "keep"),
  guides = "collect",
  highlight = NULL,
  highlight_args = list()
)
```

## Arguments

- ...:

  Tracks to stack, top to bottom, as ggplots or a single list of
  ggplots. Tracks need no coord of their own.

- windows:

  Windows to show on every track, as for
  [`coord_loci()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_loci.md),
  or `"genome"` for
  [`coord_genome()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_genome.md).
  Default: the windows of the first track that has a
  [`coord_loci()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_loci.md)
  or
  [`coord_locus()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_locus.md).

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

- heights:

  Relative heights of the tracks, default: equal

- upper_x_axis:

  Bottom x-axis style of the tracks above the bottom track: `"ticks"`
  (default), `"blank"`, or `"keep"` to leave each track's own theme.

- guides:

  Passed to
  [`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html),
  default: `"collect"` to gather legends on the right.

- highlight:

  Regions to highlight across the stack, as for
  [`highlight_stack()`](https://pblaney.github.io/ggGenomicTracks/reference/highlight_stack.md).
  Default: none.

- highlight_args:

  List of other arguments passed to
  [`highlight_stack()`](https://pblaney.github.io/ggGenomicTracks/reference/highlight_stack.md),
  e.g. `list(fill = c("#9BD08A", "#DC766F"), legend_title = "Windows")`.

## Value

A patchwork, or with `highlight`, a highlighted stack from
[`highlight_stack()`](https://pblaney.github.io/ggGenomicTracks/reference/highlight_stack.md)

## Examples

``` r
if (FALSE) { # \dontrun{
windows <- c(MYC = "chr8:126500001-129500000", IGH = "chr14:105550001-106100000")

stack_tracks(
  ggplot2::ggplot(junctions) + geom_junction(ggplot2::aes(colour = type)) + theme_genomic(),
  ggplot2::ggplot(cn) + geom_copy_number(cn_cap = 8) + theme_genomic(),
  ggplot2::ggplot() + geom_gene_track(genes_of_interest = c("MYC", "IGHM")) +
    theme_genomic(y_axis = "blank"),
  windows = windows, heights = c(1, 1, 0.5))
} # }
```
