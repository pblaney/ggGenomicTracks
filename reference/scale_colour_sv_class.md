# Junction orientation color scale

Colors for junction orientation classes from
[`annotate_junctions()`](https://pblaney.github.io/ggGenomicTracks/reference/annotate_junctions.md).
The `"jabba"` palette follows the JaBbA paper, with DEL-like junctions
in cyan and DUP-like in magenta; the `"pcawg"` palette follows the PCAWG
structural variant paper, with deletion-type in blue, duplication-type
in red, and head-to-head and tail-to-tail inversions in green and
purple.

## Usage

``` r
scale_colour_sv_class(palette = c("jabba", "pcawg"), name = "Junction", ...)
```

## Arguments

- palette:

  `"jabba"` (default) or `"pcawg"`

- name:

  Legend title, default: "Junction"

- ...:

  Other arguments passed to
  [`ggplot2::scale_colour_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html).

## Value

ggplot2 scale
