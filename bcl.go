package codeloops

import (
	"fmt"
	"strconv"
)

const (
	Pos uint = iota
	Neg
)

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
// code.
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

// CL is a structure which contains information about the ambient loop space
// such as the basis.
type CL struct {
	basis    []uint
	basisLen uint // elems in the basis
	size     uint // elems in the loop
}

// CLElem is one element of a loop. It consists of a signed vector in the
// generator space. For example, the [8,4] Hamming code gives us 16 possible
// code words, using a 4 element basis, for a total of 32 loop elements.
type CLElem struct {
	sgn uint
	vec uint
	str string
}

// NewCL returns a new code loop from a basis for a doubly even binary code.
// The basis is NOT CHECKED, because the verification is performance-heavy. To
// verify several things about the code loop, use Verify().
func NewCL(basis []uint) (cl *CL, e error) {
	cl = new(CL)
	cl.basis = basis
	cl.basisLen = uint(len(basis))
	cl.size = 1 << cl.basisLen
	return
}

// NewCLElem creates a signed entry in the loop represented by CL.
func (cl *CL) NewCLElem(x, sgn uint) (cle *CLElem, e error) {
	if x > cl.size-1 {
		e = fmt.Errorf("Bad size for element, can't make %d from %d element basis.", x, cl.basisLen)
		return
	}
	cle = new(CLElem)
	if sgn > 1 {
		e = fmt.Errorf("Bad sign value %d. (Suggestion: use the constants cl.Pos and cl.Neg)", sgn)
		return
	}
	cle.sgn = sgn
	// for each bit set in x, xor in the corresponding basis vector in the
	// code space. This embeds the entry into the larger ambient space. Eg the
	// Hamming basis allows a 4 bit x, and produces a vec in 8 bit space. Note
	// that all of this is being trivially embedded again into 64 bits so that
	// we can use simple machine arithmetic. That doesn't effect the maths.
	for idx := uint(0); idx < cl.basisLen; idx++ {
		if x&1 == 1 {
			cle.vec ^= cl.basis[idx]
		}
		x >>= 1
	}
	return
}

func (cl *CL) Size() int {
	return int(cl.size)
}

func (cl *CL) LoopElems() (cles []CLElem) {
	// We put all the positive elements first. That way, if we only want the
	// loop vectors (and not the signs) we can just pull out the first half of
	// them elems and ignore cle.sgn
	for i := uint(0); i < cl.size; i++ {
		ipos, _ := cl.NewCLElem(i, Pos)
		cles = append(cles, *ipos)
	}
	for i := uint(0); i < cl.size; i++ {
		ineg, _ := cl.NewCLElem(i, Neg)
		cles = append(cles, *ineg)
	}
	return
}

// Verify checks the supplied basis to ensure that it is a doubly even binary
// code.
func (cl *CL) Verify() (e error) {
	// TODO
	return
}

// Mul performs "multiplication" (the loop action) in the loop.
func (cle1 *CLElem) Mul(cle2 *CLElem) (res *CLElem) {
	res = new(CLElem)
	// this is the sigma function that maps from C^2 -> {0,1}
	res.sgn = cle1.sgn ^ cle2.sgn ^ (((cle1.vec & cle2.vec) >> 1) & 1)
	// this is addition in the ambient vector field. Because it has
	// characteristic 2, addition is xor. The fact that the loop vectors are
	// closed under ^ is a property of the codes. EG the 16 Hamming codewords
	// are floating in 8 bit space, but any combination of them under XOR
	// produces another of the same 16 codewords.
	res.vec = cle1.vec ^ cle2.vec
	return
}

// Mul performs "multiplication" (the loop action) in the loop.
func (cle1 *CLElem) Mul2(cle2, res *CLElem) {
	// this is the sigma function that maps from C^2 -> {0,1}
	res.sgn = cle1.sgn ^ cle2.sgn ^ (((cle1.vec & cle2.vec) >> 1) & 1)
	// this is addition in the ambient vector field. Because it has
	// characteristic 2, addition is xor. The fact that the loop vectors are
	// closed under ^ is a property of the codes. EG the 16 Hamming codewords
	// are floating in 8 bit space, but any combination of them under XOR
	// produces another of the same 16 codewords.
	res.vec = cle1.vec ^ cle2.vec
	return
}

func (c *CLElem) String() string {
	if c.str != "" {
		return c.str
	}
	sgn := ""
	if c.sgn == Neg {
		sgn = "-"
	} else {
		sgn = "+"
	}
	c.str = sgn + strconv.Itoa(int(c.vec))
	return c.str
}
