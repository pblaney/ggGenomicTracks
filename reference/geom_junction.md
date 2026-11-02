# Structural variant junction track

Draws structural variant junctions as arcs between their two
breakpoints. Pair with
[`coord_loci()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_loci.md)
to connect breakpoints in different windows, including trans junctions
between chromosomes, or
[`coord_genome()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_genome.md)
for a genome-wide view. A junction with one breakpoint out of view is
drawn as a stub at the breakpoint in view, labeled with the chromosome
of its partner.

`x` and `xend` default to the `start1` and `start2` columns of BEDPE
data, as from `devgru::read_bedpe_file(output_type = "dt")`. On its own,
the track draws arcs up from the bottom of the panel, or down from the
top, taller for junctions spanning more of the axis. To draw junctions
on a copy number track, map `y` and `yend` to the copy number at each
breakpoint, e.g. from
[`junction_copy_number()`](https://pblaney.github.io/ggGenomicTracks/reference/junction_copy_number.md),
and arcs connect those heights.

## Usage

``` r
geom_junction(
  mapping = NULL,
  data = NULL,
  arc_height = NULL,
  arc_direction = c("up", "down"),
  scale_by_span = TRUE,
  arc_shape = c("arch", "strand"),
  partner_label_size = 2.5,
  breakpoint_marks = FALSE,
  ...,
  na.rm = FALSE,
  show.legend = NA,
  inherit.aes = TRUE
)
```

## Arguments

- mapping:

  Set of aesthetic mappings created by
  [`ggplot2::aes()`](https://ggplot2.tidyverse.org/reference/aes.html),
  e.g. `colour = type`, or `y` and `yend` to anchor arcs at copy number.

- data:

  Data frame of junctions with `seqnames1`, `start1`, `seqnames2`, and
  `start2` columns.

- arc_height:

  Height of the tallest arcs as a fraction of the panel height. Default:
  0.9 on their own, 0.25 when `y` is mapped.

- arc_direction:

  Draw arcs `"up"` (default) or `"down"` from their breakpoints.

- scale_by_span:

  Make arcs spanning more of the axis taller, default: TRUE

- arc_shape:

  Shape of the arcs. `"arch"` (default) draws symmetric arches;
  `"strand"` draws arcs leaving each breakpoint toward the side it
  joins, from the `strand1` and `strand2` columns, as in gTrack:
  deletion-like junctions arch inward, duplication-like junctions loop
  outward, and inversions hook to one side, so short fold-back
  inversions show as hooks. Stubs lean the same way.

- partner_label_size:

  Size of the partner chromosome labels on stubs in mm, or 0 for none.
  Default: 2.5

- breakpoint_marks:

  Mark each breakpoint with a short tick pointing to the side it joins,
  from the `strand1` and `strand2` columns, as in PCAWG figures.
  Default: FALSE

- ...:

  Other arguments passed to
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html),
  e.g. `colour`, `linewidth`, or `alpha`.

- na.rm, show.legend, inherit.aes:

  See
  [`ggplot2::layer()`](https://ggplot2.tidyverse.org/reference/layer.html).

## Value

List of a ggplot2 layer and, unless `y` is mapped, a y-axis scale and
title, added to a plot with `+`

## Examples

``` r
if (FALSE) { # \dontrun{
junctions <- devgru::read_bedpe_file("sample.jabba.junctions.bedpe", output_type = "dt")
windows <- c(MYC = "chr8:126000001-130000000", IGH = "chr14:105000001-107000000")

ggplot2::ggplot(junctions) +
  geom_junction(ggplot2::aes(colour = type)) +
  scale_colour_sv_type() +
  coord_loci(windows) +
  theme_genomic()

# Junctions joining copy number segments
ggplot2::ggplot() +
  geom_copy_number(data = cn, cn_cap = 8) +
  geom_junction(data = junction_copy_number(junctions, cn, cn_cap = 8),
                ggplot2::aes(y = cn1, yend = cn2, colour = type)) +
  scale_colour_sv_type() +
  coord_loci(windows) +
  theme_genomic()
} # }
```
