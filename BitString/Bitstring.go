package Bitstring

import (
	// need fmt.Print for PrintBits
	"fmt"
	// need math.Ceil in NewBitstring, math.Floor in getArrayPos
	"math"
)

type Bitstring struct {
	// numBits is the number of number of bits the BitString holds
	numBits int
	// Bits is the array holding the bits. It is of length ceil(numBits/8)
	Bits []byte
}

// NewBitstring returns a bit string of n bits.
// We make ceil(numBits/8) bytes of data
func NewBitstring(numBits int) *Bitstring {
	numBytes := int(math.Ceil(float64(numBits) / float64(8.0)))
	return &Bitstring{numBits: numBits, Bits: make([]byte, numBytes)}
}

// SetBit sets the nth bit of b
func (b *Bitstring) SetBit(n int) {
	bytenum, bitnum := b.getArrayPos(n)
	b.Bits[bytenum] |= 1 << uint(bitnum)
}

// UnsetBit unsets the nth bit of b
func (b *Bitstring) UnsetBit(n int) {
	bytenum, bitnum := b.getArrayPos(n)
	b.Bits[bytenum] &= ^(1 << uint(bitnum))
}

// GetBit returns the nth bit of b
func (b *Bitstring) GetBit(n int) byte {
	bytenum, bitnum := b.getArrayPos(n)
	return (b.Bits[bytenum] >> uint(bitnum)) & 1
}

// getArrayPos returns the byte that the nth bit fits into
// as well as the bit position relative to be beginning
// of the bit string.
// For instance, if b.numBits > 13, getArrayPos(b, 5) = 0,5
//				    getArrayPos(b, 10) = 1,2.
func (b *Bitstring) getArrayPos(n int) (int, int) {
	bytenum := int(math.Floor(float64(n) / float64(8.0)))
	bitnum := n % 8
	return bytenum, bitnum
}

// PrintBits prints the bits of b
func (b *Bitstring) PrintBits() {
	for i := 0; i < int(b.numBits); i++ {
		fmt.Print(b.GetBit(i))
	}
	fmt.Println()
}
