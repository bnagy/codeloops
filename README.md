# codeloops

## About

`codeloops` is initial maths research code to implement code loops as extensions of doubly even binary codes.  Initial goal is to build the loop over the Hamming[8,4] code and then Parker loop over the Golay code. This is undergraduate code with no practical relevance for anyone except me.

## Installation

You should follow the [instructions](https://golang.org/doc/install) to
install Go, if you haven't already done so. Then:
```bash
$ go get github.com/bnagy/codeloops
```

## Usage

So far there are only tests. Run them:

```
go test
```

## License & Acknowledgements

- FOR THE CODE: BSD Style. See [LICENSE](LICENSE.md) file for details.
- FOR THE PAPER: CC-BY, see [`paper/README.md`](paper/README.md)
- FOR THE VENDORED [BITSTRING](Bitstring/README.md) CODE: _NO_ license (original code not licensed)

## TODO

Re-implement BitString, since it is fairly simple stuff and carries no license. :/

## Contributing

Fork and PR, but please don't.
