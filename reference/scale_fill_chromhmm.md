# Chromatin state fill scale

Fill colors for the 12-state ChromHMM model used by
[`geom_epilogo_col()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_epilogo_col.md).
States are matched by name with spaces, as made by
[`epilogos_long()`](https://pblaney.github.io/ggGenomicTracks/reference/epilogos_long.md),
and listed in the legend from promoters to heterochromatin.

## Usage

``` r
scale_fill_chromhmm(
  palette = c("colorblind", "epilogos"),
  name = "ChromHMM\nStates",
  ...
)
```

## Arguments

- palette:

  `"colorblind"` (default) for a color-blind friendly palette, or
  `"epilogos"` for the standard epilogos state colors.

- name:

  Legend title, default: "ChromHMM\nStates"

- ...:

  Other arguments passed to
  [`ggplot2::scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html).

## Value

ggplot2 scale
