# Bitstring

[<img src = "http://img.shields.io/badge/docs-GoDoc-red.svg">](https://godoc.org/github.com/nigelredding/Bitstring)


## Summary

This package implements the Bitstring data type. This was primarily implemented for use in github.com/nigelredding/bloom_and_cuckoo.

## Installation

To install, run

```bash
go get -u github.com/nigelredding/Bitstring
```

## Usage

To use the library, import it with

```go
import "github.com/nigelredding/Bitstring"
```

Here is an example which uses all of the exported methods in BitString:

```go
package main

import (
	"github.com/nigelredding/Bitstring"
)

func main() {
	// creates a new BitString of length 10. Uses up ceil(10/8)=2 bytes
	b := Bitstring.NewBitstring(10)

	// Prints out 0000000000, as all bits are initially 0
	b.PrintBits()

	// set all the bits to 1
	for i := 0; i < 10; i++ {
		b.SetBit(i)
	}

	// Prints out 1111111111, as expected
	b.PrintBits()

	// set the first 5 bits to 0
	for i := 0; i < 5; i++ {
		b.UnsetBit(i)
	}

	// Prints out 0000011111, as expected
	b.PrintBits()

	// set the last bit to 0
	b.UnsetBit(9)

	// Prints out 0000011110, as expected
	b.PrintBits()
}
```
