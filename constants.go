package codeloops

const (
	Pos uint = iota
	Neg
)

const RandomTheta = 0xffffffff

// HammingBasis is a basis for the [8,4] Hamming code, which generates the 16
// possible codewords in that space.
var HammingBasis = []uint{
	0x87, // 1000 | 0111
	0x4b, // 0100 | 1011
	0x2d, // 0010 | 1101
	0x1e, // 0001 | 1110
}

// 100000000000100111110001 = 0x8009f1
// 010000000000010011111010 = 0x4004fa
// 001000000000001001111101 = 0x20027d
// 000100000000100100111110 = 0x10093e
// 000010000000110010011101 = 0x80c9d
// 000001000000111001001110 = 0x40e4e
// 000000100000111100100101 = 0x20f25
// 000000010000111110010010 = 0x10f92
// 000000001000011111001001 = 0x87c9
// 000000000100001111100110 = 0x43e6
// 000000000010010101010111 = 0x2557
// 000000000001101010101011 = 0x1aab

// GolayBasis is a basis for the binary Golay code, a [24,12] error correcting
// code. This is a 'standard' basis which appears in multiple sources.
var GolayBasis = []uint{
	0x8009f1,
	0x4004fa,
	0x20027d,
	0x10093e,
	0x80c9d,
	0x40e4e,
	0x20f25,
	0x10f92,
	0x87c9,
	0x43e6,
	0x2557,
	0x1aab,
}

// GolaySplitBasis is a different basis which has the property that the first
// 6 vectors split into a group, the subsequent 5 also split into a group, and
// those two subspaces meet only in the zero vector. The final vector is
// arbitrary to increase the span to 4096.
var GolaySplitBasis = []uint{
	0x8009f1,
	0x4004fa,
	0x10093e,
	0x80c9d,
	0x2557,
	0x27d200,
	0x20027d,
	0x40e4e,
	0x810663,
	0x808e38,
	0x804a17,
	0x6c1eff, // 0xb1f6f3,
}
