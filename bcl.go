package codeloops

import (
	"fmt"
	"runtime"
	"strconv"
	"strings"
	"sync"
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

func (cl *CL) PrintLoopElems() {
	fmt.Print("----------\n")
	for i, cle := range cl.LoopElems() {
		fmt.Printf("%d: %s\n", i, &cle)
	}
	fmt.Print("----------\n")
}

// VerifyBasis checks the supplied basis to ensure that it is a doubly even binary
// code.
func (cl *CL) VerifyBasis() (e error) {
	for _, cle := range cl.LoopElems() {
		if bitWeight(cle.vec)%4 != 0 {
			e = fmt.Errorf("Bad vector %b, bitweight not a multiple of 4.", cle.vec)
			return
		}
	}
	return
}

func (cl *CL) VerifyMoufang() (e error) {
	e = cl.checkMoufangParallel()
	return
}

func (cl *CL) VerifyMoufangSlow() (e error) {
	e = cl.checkMoufangParallelSlow()
	return
}

func sigma(x, y uint) uint {
	return (bitWeight(x&y) >> 1) & 1
}

// Mul performs "multiplication" (the loop action) in the loop.
func (res *CLElem) Mul(x, y *CLElem) *CLElem {
	// this is the sigma function that maps from C^2 -> {0,1}
	res.sgn = x.sgn ^ y.sgn ^ sigma(x.vec, y.vec)
	// this is addition in the ambient vector field. Because it has
	// characteristic 2, addition is xor. The fact that the loop vectors are
	// closed under ^ is a property of the codes. EG the 16 Hamming codewords
	// are floating in 8 bit space, but any combination of them under XOR
	// produces another of the same 16 codewords.
	res.vec = x.vec ^ y.vec
	return res
}

func (c *CLElem) String() string {
	sgn := ""
	if c.sgn == Neg {
		sgn = "-"
	} else {
		sgn = "+"
	}
	c.str = sgn + strconv.Itoa(int(c.vec))
	return c.str
}

type WorkUnit struct {
	BaseCL   *CL
	Elems    []CLElem // readonly, should be safe
	StartIdx uint
	StopIdx  uint
	Failures chan<- []*CLElem
	Wg       *sync.WaitGroup
}

func moufangParallelWorker(work *WorkUnit) {

	defer work.Wg.Done()

	// Preallocate the memory for the multiplication elems
	x, y, z := new(CLElem), new(CLElem), new(CLElem)

	// Iterate! Treat a bit string as 3 concatenated element indices, and then
	// use those to select the appropriate real element out of the work.Elems
	// slice.

	// Create masks to pull the element indicies from the counter
	elemSz := int(work.BaseCL.basisLen)
	zeroChunk := strings.Repeat("0", elemSz)
	oneChunk := strings.Repeat("1", elemSz)
	yMaskStr := zeroChunk + oneChunk + zeroChunk
	zMaskStr := zeroChunk + zeroChunk + oneChunk
	yMask, _ := strconv.ParseUint(yMaskStr, 2, 0)
	zMask, _ := strconv.ParseUint(zMaskStr, 2, 0)

	for combinedElems := work.StartIdx; combinedElems < work.StopIdx; combinedElems++ {
		x = &work.Elems[combinedElems>>uint(elemSz*2)]             // top elemSz bits as an index
		y = &work.Elems[(combinedElems&uint(yMask))>>uint(elemSz)] // middle elemSz bits ...
		z = &work.Elems[combinedElems&uint(zMask)]                 // etc

		// This is a Moufang identity expressed in XOR. If:
		// sigma(x,z)sigma(y,xz)sigma(x,xyz) == sigma(x,y)sigma(xy,x)sigma(y,z)
		// then the Moufang property holds.
		if sigma(x.vec, z.vec)^
			sigma(y.vec, x.vec^z.vec)^
			sigma(x.vec, x.vec^y.vec^z.vec)^
			sigma(x.vec, y.vec)^
			sigma(x.vec^y.vec, x.vec)^
			sigma(y.vec, z.vec) != 0 {
			work.Failures <- []*CLElem{x, y, z}
			return // Assuming all the other workers also exit early, we'll be OK.
			// Worst case, we run through ~all the tests, which is still twice
			// as fast as using an abort channel via a select() loop here.
		}
	}
}

func (cl *CL) checkMoufangParallel() (e error) {

	elems := cl.LoopElems()

	wg := &sync.WaitGroup{}
	failures := make(chan []*CLElem)
	done := make(chan struct{})

	go func() {
		// This goroutine will wait for all the workers to finish, then close
		// the done channel, which aborts the main select loop.
		wg.Wait()
		close(done)
	}()

	// create and dispatch the work units. Each worker gets assigned a chunk
	// of the whole range divided by the number of logical CPUs.
	elemSz := cl.basisLen
	lim := 1 << (elemSz * 3)
	step := lim / runtime.NumCPU() // more or less
	for cpu := 0; cpu < runtime.NumCPU(); cpu++ {
		// create this work unit, and start the worker
		wg.Add(1)
		thisStart := uint(step * cpu)
		thisStop := uint(0)
		if cpu+1 == runtime.NumCPU() { // last cpu gets the slack
			thisStop = uint(lim)
		} else {
			thisStop = uint((cpu + 1) * step)
		}
		wu := &WorkUnit{cl, elems, thisStart, thisStop, failures, wg}
		go moufangParallelWorker(wu)
	}

loop:
	for {
		select {
		case fail := <-failures:
			// there could be many errors, we just keep the last one.
			e = fmt.Errorf("Failed Moufang Identity for %s %s %s", fail[0].String(), fail[1].String(), fail[2].String())
		case <-done:
			// will fire when the done chan is closed (all workers have returned)
			break loop
		}
	}
	return
}

// Older and slower, but uses full Mul() calls. Useful in case of bugs.
func moufangParallelWorkerSlow(work *WorkUnit) {

	defer work.Wg.Done()

	// Preallocate the memory for the multiplication elems
	// Preallocate the memory for the multiplication elems
	x, y, z := new(CLElem), new(CLElem), new(CLElem)
	zy, x_zy, z__x_zy := new(CLElem), new(CLElem), new(CLElem)
	zx, zx_z, zx_z__y := new(CLElem), new(CLElem), new(CLElem)

	// Iterate! Treat a bit string as 3 concatenated element indices, and then
	// use those to select the appropriate real element out of the work.Elems
	// slice.

	// Create masks to pull the element indicies from the counter
	elemSz := work.BaseCL.basisLen + 1
	zeroChunk := strings.Repeat("0", int(elemSz))
	oneChunk := strings.Repeat("1", int(elemSz))
	yMaskStr := zeroChunk + oneChunk + zeroChunk
	zMaskStr := zeroChunk + zeroChunk + oneChunk
	yMask, _ := strconv.ParseUint(yMaskStr, 2, 0)
	zMask, _ := strconv.ParseUint(zMaskStr, 2, 0)

	for combinedElems := work.StartIdx; combinedElems < work.StopIdx; combinedElems++ {
		// verify Moufang identity:
		// z(x(zy)) = ((zx)z)y
		x = &work.Elems[combinedElems>>(elemSz*2)]           // top 5 bits as an index
		y = &work.Elems[(combinedElems&uint(yMask))>>elemSz] // middle 5 bits ...
		z = &work.Elems[combinedElems&uint(zMask)]           // etc

		// LHS 1
		zy.Mul(z, y)
		x_zy.Mul(x, zy)
		z__x_zy.Mul(z, x_zy)

		//RHS 1
		zx.Mul(z, x)
		zx_z.Mul(zx, z)
		zx_z__y.Mul(zx_z, y)

		// The vector ops are associative, the only thing that can mismatch is the sign.
		if z__x_zy.sgn != zx_z__y.sgn {
			work.Failures <- []*CLElem{x, y, z}
		}

	}

}

func (cl *CL) checkMoufangParallelSlow() (e error) {

	elems := cl.LoopElems()

	wg := &sync.WaitGroup{}
	failures := make(chan []*CLElem)
	done := make(chan struct{})

	go func() {
		// This goroutine will wait for all the workers to finish, then close
		// the done channel, which aborts the main select loop.
		wg.Wait()
		close(done)
	}()

	// create and dispatch the work units. Each worker gets assigned a chunk
	// of the whole range divided by the number of logical CPUs.
	elemSz := cl.basisLen + 1
	lim := 1 << (elemSz * 3)
	step := lim / runtime.NumCPU() // more or less
	for cpu := 0; cpu < runtime.NumCPU(); cpu++ {
		// create this work unit, and start the worker
		wg.Add(1)
		thisStart := uint(step * cpu)
		thisStop := uint(0)
		if cpu+1 == runtime.NumCPU() { // last cpu gets the slack
			thisStop = uint(lim)
		} else {
			thisStop = uint((cpu + 1) * step)
		}
		wu := &WorkUnit{cl, elems, thisStart, thisStop, failures, wg}
		go moufangParallelWorkerSlow(wu)
	}

	// Wait for the work to finish, or bail immediately if one of the workers
	// reports a failed triple.
loop:
	for {
		select {
		case fail := <-failures:
			// there could be many errors, we just keep the last one.
			e = fmt.Errorf("Failed Moufang Identity for %s %s %s", fail[0].String(), fail[1].String(), fail[2].String())
		case <-done:
			// will fire when the done chan is closed (all workers have returned)
			break loop
		}
	}
	return
}
