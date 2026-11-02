# Build an integer copy number scale

The palette depends on the copy numbers present, so it is built from the
scale limits when values are mapped rather than from the number of
levels. For internal use only.

## Usage

``` r
copynumber_scale(aesthetics, name, ..., na.value = "black")
```

## Arguments

- aesthetics:

  Aesthetic of the scale

- name:

  Legend title, default: "Copy\nNumber"

- ...:

  Other arguments passed to
  [`ggplot2::discrete_scale()`](https://ggplot2.tidyverse.org/reference/discrete_scale.html),
  e.g. `na.value`.
