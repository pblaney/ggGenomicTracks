# Position gene labels within the panel at draw time

Places each label at its requested gene end if it fits in the panel,
otherwise at the opposite gene end, otherwise centered above the visible
part of the gene. Called by the label grob once the panel size is known.
For internal use only.

## Usage

``` r
place_gene_labels(label_grob)
```

## Arguments

- label_grob:

  gTree of class `gene_label_grob` built by `GeomGeneLabel`
