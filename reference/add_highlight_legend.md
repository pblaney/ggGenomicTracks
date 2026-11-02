# Add a legend of highlighted regions to a stack

Builds the legend from a throwaway plot of the regions and adds it below
the collected legends of the stack, or in a new column to the right if
legends are not collected. For internal use only.

## Usage

``` r
add_highlight_legend(stack_gtable, regions, alpha, colour, legend_title)
```

## Arguments

- stack_gtable:

  Gtable of the stack

- regions:

  Standardized regions from
  [`as_highlight_regions()`](https://pblaney.github.io/ggGenomicTracks/reference/as_highlight_regions.md)

- alpha:

  Opacity of the bands, default: 0.3

- colour:

  Outline color of the bands, default: none

- legend_title:

  Title of the legend, default: "Highlight"
