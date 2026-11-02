# Lay out genomic windows on a linear axis

Assigns each window a span of the linear axis. The first window starts
at its own genomic start, so a single window keeps genomic coordinates
unchanged. For internal use only.

## Usage

``` r
locus_layout(windows, widths = "proportional", gap = 0.02)
```

## Arguments

- windows:

  Data frame of windows from
  [`as_locus_windows()`](https://pblaney.github.io/ggGenomicTracks/reference/as_locus_windows.md)

- widths:

  How the axis is shared between windows: `"proportional"` (default) to
  each window's length in bp, or `"equal"` so small windows stay
  readable next to large ones, each then at its own scale.

- gap:

  Space between windows as a fraction of the total axis, default: 0.02

## Value

Data frame of windows with `chrom`, `chrom_key` (no `chr` prefix),
`start`, `end`, `label`, `scale` (axis units per bp), `lin_start`, and
`lin_end`
