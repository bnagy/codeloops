# codeloops

## About

`codeloops` is initial maths research code to implement code loops as extensions of doubly even binary codes. For further details, consult the [paper](paper/README.md).

## Installation

You should follow the [instructions](https://golang.org/doc/install) to
install Go, if you haven't already done so. Then:
```bash
$ go get github.com/bnagy/codeloops
```

## Usage

Run tests:

```
go test
```

For descriptions of the various utility programs (some of which are really not for external consumption), see [`cmd/README.md`](cmd/README.md).

## License & Acknowledgements

- FOR THE CODE: BSD Style. See [LICENSE](LICENSE.md) file for details.
- FOR THE PAPER: CC-BY, see [`paper/README.md`](paper/README.md)
- FOR THE VENDORED [BITSTRING](BitString/README.md) CODE: _NO_ license (original code not licensed)

## TODO

Re-implement BitString, since it is fairly simple stuff and carries no license. :/

## Contributing

Fork and PR. Contact paper authors for other queries.
