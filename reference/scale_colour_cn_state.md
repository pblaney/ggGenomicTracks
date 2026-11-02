# Copy number state color scales

Colors for the JaBbA copy number states, from homozygous deletion in
dark blue to amplification in dark red, with neutral in gray. Use the
fill scale with
[`geom_copy_number()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_copy_number.md)
or with
[`geom_bed_track()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_bed_track.md)
to show states across many samples.

## Usage

``` r
scale_colour_cn_state(name = "CN State", ...)

scale_fill_cn_state(name = "CN State", ...)
```

## Arguments

- name:

  Legend title, default: "CN State"

- ...:

  Other arguments passed to
  [`ggplot2::scale_colour_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html)
  or
  [`ggplot2::scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html),
  e.g. `values` to change colors.

## Value

ggplot2 scale
