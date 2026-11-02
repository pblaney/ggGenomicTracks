# Classify junctions by orientation and complex event

Adds two columns for coloring
[`geom_junction()`](https://pblaney.github.io/ggGenomicTracks/reference/geom_junction.md):

- `sv_class`: orientation from the breakpoint strands, `"DEL-like"`
  (+/-), `"DUP-like"` (-/+), `"INV-h2h"` (+/+), or `"INV-t2t"` (-/-),
  for both cis and trans junctions. Use with
  [`scale_colour_sv_class()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_colour_sv_class.md).

- `event`: the JaBbA complex event the junction belongs to, from the
  event columns (`chromothripsis`, `chromoplexy`, `templated_insertion`,
  `rigma`, `pyrgo`, `tyfonas`, `ecdna`, `bfb`), or `"Simple"`; and
  `event_id` with the event label, e.g. `"CPX-4"`. Use with
  [`scale_colour_jabba_event()`](https://pblaney.github.io/ggGenomicTracks/reference/scale_colour_jabba_event.md).

## Usage

``` r
annotate_junctions(junctions)
```

## Arguments

- junctions:

  Data frame of junctions with `strand1` and `strand2` columns and
  optionally JaBbA event columns, e.g. from
  `devgru::read_bedpe_file(output_type = "dt")` of a JaBbA
  `junctions.bedpe`.

## Value

`junctions` as a data frame with `sv_class`, `event`, and `event_id`
added
