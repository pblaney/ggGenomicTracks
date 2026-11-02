# Locus coordinate and scale bar annotation

Adds a scale bar and the window coordinates to a track panel. The bar
length is chosen from the visible x-range (a round 1, 2, or 5 x 10^n bp
up to 20% of the width), so it adapts to any zoom. With
[`coord_locus()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_locus.md)
the coordinates show the exact window, e.g.
"chr3:187,729,001-187,756,500 (27,500 bp)"; otherwise the visible range
is shown without a chromosome. With several windows from
[`coord_loci()`](https://pblaney.github.io/ggGenomicTracks/reference/coord_loci.md),
only the scale bar is drawn, and nothing is drawn when the windows have
equal widths, as each then has its own scale.

Add to any track, or recreate a standalone scale track with
`ggplot() + annotation_locus_scale(location = "middle") + coord_locus(win) + theme_void()`.

## Usage

``` r
annotation_locus_scale(
  location = c("top", "middle", "bottom"),
  bar_side = c("right", "left"),
  show_coordinates = TRUE,
  colour = "black",
  linewidth = 0.8,
  text_size = 3.4,
  pad = 0.02
)
```

## Arguments

- location:

  Vertical placement in the panel, one of `"top"`, `"middle"`, or
  `"bottom"`.

- bar_side:

  Side of the panel for the scale bar, `"right"` or `"left"`. The
  coordinates are placed on the opposite side.

- show_coordinates:

  Show the window coordinates, default: TRUE

- colour:

  Color of the bar and text.

- linewidth:

  Width of the bar in mm, default: 0.8

- text_size:

  Size of the text in mm, default: 3.4

- pad:

  Distance from the panel edges as a fraction of panel width, default:
  0.02

## Value

ggplot2 layer

## Examples

``` r
if (FALSE) { # \dontrun{
ggplot2::ggplot() +
  geom_bed_track(data = chrom_states, ggplot2::aes(y = sample, fill = segment_label)) +
  annotation_locus_scale() +
  coord_locus("chr3:187734001-187751500", flank = 5000) +
  theme_genomic()
} # }
```
