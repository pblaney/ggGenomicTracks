# Structural variant type color scale

Colors for deletions (`DEL`), duplications (`DUP`), inversions (`INV`),
and translocations (`TRA`).

## Usage

``` r
scale_colour_sv_type(name = "SV Type", ...)
```

## Arguments

- name:

  Legend title, default: "SV Type"

- ...:

  Other arguments passed to
  [`ggplot2::scale_colour_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html),
  e.g. `values` to change colors.

## Value

ggplot2 scale
