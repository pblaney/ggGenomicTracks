# Blank the lowest axis label

Wraps a label function so the lowest visible break is left unlabeled
while its tick is kept. Useful in tightly stacked tracks, where the
lowest y-axis label of one track collides with the highest label of the
track below it.

## Usage

``` r
label_drop_first(labels = scales::label_number())
```

## Arguments

- labels:

  Label function to wrap, default:
  [`scales::label_number()`](https://scales.r-lib.org/reference/label_number.html)

## Value

Label function for the `labels` argument of a scale

## Examples

``` r
if (FALSE) { # \dontrun{
ggplot2::ggplot(gistic) +
  geom_ridge_track(ggplot2::aes(y = G.score)) +
  ggplot2::scale_y_continuous(labels = label_drop_first()) +
  coord_locus("chr13:48303751-48481890", flank = 1e6) +
  theme_genomic(x_axis = "ticks")
} # }
```
