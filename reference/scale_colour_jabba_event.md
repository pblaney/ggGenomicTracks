# JaBbA complex event color scale

Colors for the complex event classes from
[`annotate_junctions()`](https://pblaney.github.io/ggGenomicTracks/reference/annotate_junctions.md),
with junctions outside complex events in light gray so events stand out.

## Usage

``` r
scale_colour_jabba_event(name = "Event", ...)
```

## Arguments

- name:

  Legend title, default: "Event"

- ...:

  Other arguments passed to
  [`ggplot2::scale_colour_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html).

## Value

ggplot2 scale
