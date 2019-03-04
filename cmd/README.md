# codeloops/cmd

## About

This project did two main things:
- Provided a concrete reference implementation of the Griess 'algorithm' (or proof) from [Gri86]
- Showed that it is possible to perform multiplication in the loop using a preconstructed set of 'alpha' data, without the need to precalculate the entire code loop. This technique may make it easier to work with much larger loops.

The `cmd` contains three utility programs that were used during the research and preparation of the paper.

`alpha_pics` was used to create visualisations of assorted alpha squares.

`loop_pics` draws the entire Parker loop, using our special basis.

`partition` partitions the subspace by recursively choosing two bases of specified lengths and with configurable partitioning criteria. This was the tool that was used to find the 'special' bases for Golay which allowed for the construction of the alpha tables.

`subviz` visualises the subspaces formed by extending subspaces of the full code (mostly used with the Golay code). Some of these spaces split and form associative loops (ie groups), while some do not. I can't remember why I was obsessed with this.

## Contributing

Fork and PR.
