# Integer copy number color scales

Colors for integer copy number from
[`devgru::copynumber_palettier()`](https://pblaney.github.io/devgru/reference/copynumber_palettier.html):
blue for 0 and 1, gray for 2, and orange to dark red from 3 up to the
highest copy number in the plot. Map a discrete copy number, e.g.
`aes(fill = factor(total_cn))` with
[`geom_bed_track()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_bed_track.md);
[`geom_copy_number()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_copy_number.md)
uses the fill scale by default.

## Usage

``` r
scale_colour_copynumber(name = "Copy\nNumber", ...)

scale_fill_copynumber(name = "Copy\nNumber", ...)
```

## Arguments

- name:

  Legend title, default: "Copy\nNumber"

- ...:

  Other arguments passed to
  [`ggplot2::discrete_scale()`](https://ggplot2.tidyverse.org/reference/discrete_scale.html),
  e.g. `na.value`.

## Value

ggplot2 scale
