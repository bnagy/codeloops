package codeloops

import (
	"fmt"
	"runtime"
	"sort"
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
	for bitPos := uint(0); bitPos < cl.basisLen; bitPos++ {
		if x&1 == 1 {
			cle.vec ^= cl.basis[bitPos]
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
	for _, vec := range cl.VectorSpace() {
		cles = append(cles, CLElem{sgn: Pos, vec: vec})
	}
	for _, vec := range cl.VectorSpace() {
		cles = append(cles, CLElem{sgn: Neg, vec: vec})

	}
	// for i := uint(0); i < cl.size; i++ {
	// 	ipos, _ := cl.NewCLElem(i, Pos)
	// 	cles = append(cles, *ipos)
	// }
	// for i := uint(0); i < cl.size; i++ {
	// 	ineg, _ := cl.NewCLElem(i, Neg)
	// 	cles = append(cles, *ineg)
	// }
	return
}

func (cl *CL) VectorSpace() (vecs []uint) {
	x := uint(0)
	set := map[uint]struct{}{}
	for i := uint(0); i < cl.size; i++ {
		x = i
		vec := uint(0)
		for bitPos := uint(0); bitPos < cl.basisLen; bitPos++ {
			if x&1 == 1 {
				vec ^= cl.basis[bitPos]
			}
			x >>= 1
		}
		set[vec] = struct{}{}
	}
	for vec, _ := range set {
		vecs = append(vecs, vec)
	}
	sort.Slice(vecs, func(i, j int) bool { return vecs[i] < vecs[j] })
	return
}

func (cl *CL) PrintBasis() {
	fmt.Print("----------\n")
	for i, vec := range cl.basis {
		fmt.Printf("%.2d: %.8b (%d) \n", i, vec, vec)
	}
	fmt.Print("----------\n")
}

func (cl *CL) PrintLoopElems() {
	fmt.Print("----------\n")
	for i, cle := range cl.LoopElems() {
		fmt.Printf("%.2d: %.8b (%s) \n", i, cle.vec, cle.String())
	}
	fmt.Print("----------\n")
}

func (cl *CL) PrintVectorSpace() {
	fmt.Print("----------\n")
	for i, vec := range cl.VectorSpace() {
		fmt.Printf("%.2d: %.8b (%d) \n", i, vec, vec)
	}
	fmt.Print("----------\n")
}

// VerifyBasis checks the supplied basis to ensure that it is a doubly even binary
// code.
func (cl *CL) VerifyBasis() (e error) {
	for _, vec := range cl.VectorSpace() {
		if BitWeight(vec)%4 != 0 {
			e = fmt.Errorf("Bad vector %b, bitweight not a multiple of 4.", vec)
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

// TODO this is rubbish
func sigma(x, y uint) uint {
	return (BitWeight(x&y) >> 1) & 1
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

func makeJoinFunc(len uint) func(uint, uint) uint {
	f := func(a, b uint) uint {
		return (a << len) | b
	}
	return f
}

func (cl *CL) buildTheta2() (map[uint]uint, error) {

	join := makeJoinFunc(8) // HAX!
	theta := make(map[uint]uint)
	b0 := cl.basis[0]

	// basic assumptions
	theta[join(b0, b0)] = (BitWeight(b0) / 4) % 2
	theta[join(0, b0)] = 0
	theta[join(b0, 0)] = 0
	theta[0] = 0

	Vk := []uint{0, b0}

	for _, bi := range cl.basis[1:] {

		// So that Vk+1 == Vk + Wk
		Wk := []uint{}
		for _, v := range Vk {
			Wk = append(Wk, bi^v)
		}

		fmt.Printf("Looping with %x. Vk: %#v Wk: %#v\n", bi, Vk, Wk)
		// D1
		for _, v := range Vk {
			fmt.Printf("D1-1 - set %.4x to %d\n", join(bi, v), 0)
			theta[join(bi, v)] = 0
			fmt.Printf("D1-2 - set %.4x to %d\n", join(v, bi), (BitWeight(v&bi)/2)%2)
			theta[join(v, bi)] = (BitWeight(v&bi) / 2) % 2
		}

		// D2
		for _, v := range Vk {
			a, ok := theta[join(bi, v)]
			if !ok {
				return nil, fmt.Errorf("Missing entry for %x in theta in D2", join(bi, v))
			}
			fmt.Printf("D2-1 - set %.4x to %d\n", join(bi, bi^v), (BitWeight(bi)/4+a)%2)
			theta[join(bi, bi^v)] = (BitWeight(bi)/4 + a) % 2
			fmt.Printf("D2-2 - set %.4x to %d\n", join(bi^v, bi), (BitWeight(bi&(bi^v))/2+theta[join(bi, bi^v)])%2)
			theta[join(bi^v, bi)] = (BitWeight(bi&(bi^v))/2 + theta[join(bi, bi^v)]) % 2
		}

		// D3
		for _, v := range Vk {
			for _, w := range Wk {
				a, ok := theta[join(v, bi)]
				if !ok {
					return nil, fmt.Errorf("Missing entry 1 for %x in theta in D3", join(v, bi))
				}
				b, ok := theta[join(v, bi^w)]
				if !ok {
					return nil, fmt.Errorf("Missing entry 2 for %x in theta in D3", join(v, bi^w))
				}
				c, ok := theta[join(w, bi)]
				if !ok {
					return nil, fmt.Errorf("Missing entry 3 for %x in theta in D3", join(w, bi))
				}
				fmt.Printf("D3- set %.4x to %d\n", join(w, bi^v), (BitWeight(v&w)/2+a+b+c)%2)

				theta[join(w, bi^v)] = (BitWeight(v&w)/2 + a + b + c) % 2
			}
		}

		// D4
		for _, w := range Wk {
			for _, v := range Vk {
				a, ok := theta[join(w, v^w)]
				if !ok {
					return nil, fmt.Errorf("Missing entry for %x in theta in D4", join(w, v^w))
				}
				fmt.Printf("w: %.4x v: %.4x w^v: %.4x a: %d (theta(%.4x) bw w/4 %d \n", w, v, w^v, a, join(w, v^w), BitWeight(w)/4)
				fmt.Printf("D4-1 set %.4x to %d\n", join(w, v), ((BitWeight(w)/4)+a)%2)
				fmt.Printf("D4-2 set %.4x to %d\n", join(v, w), (BitWeight(v&w)/2+theta[join(w, v)])%2)

				theta[join(w, v)] = ((BitWeight(w) / 4) + a) % 2
				theta[join(v, w)] = (BitWeight(v&w)/2 + theta[join(w, v)]) % 2
			}
		}

		Vk = append(Vk, Wk...)

	}

	return theta, nil
}

func (cl *CL) buildTheta() map[uint]uint {
	// Theta is a map from 2 code elems to 0 or 1 - ie if we have 16 code
	// elems there are 256 possible theta inputs. A basis of length 4
	// generates 16 vectors in the code space, so we need (2**basisLen)**2 map
	// entries.

	// this will need to be refactored into a bitstring
	join := makeJoinFunc(8) // HAX!
	theta := make(map[uint]uint)
	v0 := cl.basis[0]

	// basic assumptions
	theta[join(v0, v0)] = (BitWeight(v0) / 4) % 2
	theta[join(0, v0)] = 0
	theta[join(v0, 0)] = 0
	theta[0] = 0

	workingBasis := []uint{v0}
	tmpCL, _ := NewCL(workingBasis)
	workingSpace := tmpCL.VectorSpace()

	// Build up the theta map. We know that certain identities must hold for
	// the theta properties to be satisfied, so we solve those constraints by
	// adding basis vectors one by one and building the map inductively
	// (iteratively)
	for _, vi := range cl.basis[1:] {

		// D1
		for i := 0; i < len(workingSpace); i++ {
			x := workingSpace[i]
			theta[join(vi, x)] = 0 // an arbitrary choice. The rest must follow.
			theta[join(x, vi)] = (BitWeight(x&vi) / 2) % 2
		}

		// D2
		for i := 0; i < len(workingSpace); i++ {
			x := workingSpace[i]
			theta[join(vi, vi^x)] = ((BitWeight(vi)) / 4) % 2
			theta[join(vi^x, vi)] = (BitWeight(vi)/4 + BitWeight(vi&(vi^x))/2) % 2
		}

		// D3
		for i := 0; i < len(workingSpace); i++ {
			for j := 0; j < len(workingSpace); j++ {
				x, y := workingSpace[i], workingSpace[j]

				thetaxy, ok := theta[join(x, y)]
				if !ok {
					panic("missing entry in theta")
				}

				theta[join(vi^x, vi^y)] = (BitWeight(y&(vi^x))/2 +
					BitWeight(y&vi)/2 +
					thetaxy +
					BitWeight(vi)/4 +
					BitWeight(vi&(vi^x))/2) % 2
			}
		}

		// D4
		for i := 0; i < len(workingSpace); i++ {
			for j := 0; j < len(workingSpace); j++ {
				x, y := workingSpace[i], workingSpace[j]
				thetvixvixy, ok := theta[join(vi^x, vi^x^y)]
				if !ok {
					panic("missing entry in theta")
				}
				theta[join(vi^x, y)] = (BitWeight(vi^x)/4 + thetvixvixy) % 2
				theta[join(y, vi^x)] = (BitWeight(y&(vi^x))/2 + BitWeight(vi^x)/4 + thetvixvixy) % 2
			}
		}

		workingBasis = append(workingBasis, vi)
		tmpCL, _ = NewCL(workingBasis)
		workingSpace = tmpCL.VectorSpace()

	}
	return theta

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

	go func() {
		// This goroutine will wait for all the workers to finish, then close
		// the done channel, which aborts the main select loop.
		wg.Wait()
		close(done)
	}()

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
