# Convert epilogos scores to long format

Converts epilogos score output, with one score column per chromatin
state, to one row per bin and state with `seqnames`, `start`, `end`,
`state`, and `score` columns. Underscores in state names are replaced
with spaces and states are ordered as their columns.

## Usage

``` r
epilogos_long(data, state_names = chromhmm_12_states)
```

## Arguments

- data:

  Data frame or GRanges of epilogos scores. The first three columns of a
  data frame are taken as the chromosome, start, and end; every other
  column except `strand` and `width` is a state score.

- state_names:

  Names of the states in the order of the score columns, default: the
  12-state ChromHMM model in epilogos order. Use `NULL` to keep the
  existing column names.

## Value

data.frame in long format
