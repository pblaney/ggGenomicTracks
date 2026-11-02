# Cytoband annotations per chromosome on hg38

Annotation of cytoband name and G-banding per chromosome

## Usage

``` r
cytobands_hg38
```

## Format

A `GRanges` object with 862 ranges across chr1-22, X, and Y and 2
metadata columns:

- name:

  Cytoband name, e.g. `p36.33`

- gieStain:

  Giemsa stain result, e.g. `gneg`, `gpos25`, `acen`

## Source

UCSC Genome Browser hg38 `cytoBand` table
