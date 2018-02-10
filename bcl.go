package codeloops

import (
	"fmt"
	"github.com/nigelredding/BitString"
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
	theta    *Bitstring.Bitstring
	vs       []uint
	vm       map[uint]uint
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
	cl.theta = Bitstring.NewBitstring(int(cl.size * cl.size))
	cl.vs = cl.VectorSpace()
	cl.vm = cl.VectorIdxMap()
	return
}

// NewElemFromIdx creates a signed entry in the loop represented by CL.
func (cl *CL) NewElemFromIdx(idx, sgn uint) (cle *CLElem, e error) {
	if idx > cl.size-1 {
		e = fmt.Errorf("Bad index for vector. Can't make %d from %d element basis.", idx, cl.basisLen)
		return
	}
	cle = new(CLElem)
	if sgn > 1 {
		e = fmt.Errorf("Bad sign value %d. (Suggestion: use the constants cl.Pos and cl.Neg)", sgn)
		return
	}
	cle.sgn = sgn
	// for each bit set in idx, xor in the corresponding basis vector in the
	// code space. This embeds the entry into the larger ambient space. Eg the
	// Hamming basis allows a 4 bit x, and produces a vec in 8 bit space. Note
	// that all of this is being trivially embedded again into 64 bits so that
	// we can use simple machine arithmetic. That doesn't effect the maths.
	for bitPos := uint(0); bitPos < cl.basisLen; bitPos++ {
		if idx&1 == 1 {
			cle.vec ^= cl.basis[bitPos]
		}
		idx >>= 1
	}
	return
}

// NewElem creates a signed entry in the loop represented by CL.
func (cl *CL) NewElem(vec, sgn uint) (cle *CLElem, e error) {

	_, ok := cl.vm[vec]
	if !ok {
		e = fmt.Errorf("Vector %x is not in the underlying space", vec)
		return
	}
	cle = new(CLElem)
	if sgn > 1 {
		e = fmt.Errorf("Bad sign value %d. (Suggestion: use the constants cl.Pos and cl.Neg)", sgn)
		return
	}
	cle.vec = vec
	cle.sgn = sgn
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
		cle, _ := cl.NewElem(vec, Pos)
		cles = append(cles, *cle)
	}
	for _, vec := range cl.VectorSpace() {
		cle, _ := cl.NewElem(vec, Neg)
		cles = append(cles, *cle)
	}
	return
}

func (cl *CL) VectorSpace() (vecs []uint) {

	if cl.vs != nil {
		return cl.vs
	}

	x := uint(0)
	for i := uint(0); i < cl.size; i++ {
		x = i
		vec := uint(0)
		for bitPos := uint(0); bitPos < cl.basisLen; bitPos++ {
			if x&1 == 1 {
				vec ^= cl.basis[bitPos]
			}
			x >>= 1
		}
		vecs = append(vecs, vec)
	}
	return
}

func (cl *CL) VectorIdxMap() (m map[uint]uint) {

	if cl.vm != nil {
		return cl.vm
	}

	m = make(map[uint]uint)
	for i, v := range cl.VectorSpace() {
		m[v] = uint(i)
	}
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

// Mul performs "multiplication" (the loop action) in the loop.
func (cl *CL) Mul(x, y, res *CLElem) (*CLElem, error) {
	// this is the sigma function that maps from C^2 -> {0,1}
	if cl.theta == nil {
		cl.BuildTheta()
	}
	t, err := cl.ThetaByVec(x.vec, y.vec)
	if err != nil {
		return nil, err
	}
	res.sgn = x.sgn ^ y.sgn ^ t
	// this is addition in the ambient vector field. Because it has
	// characteristic 2, addition is xor. The fact that the loop vectors are
	// closed under ^ is a property of the codes. EG the 16 Hamming codewords
	// are floating in 8 bit space, but any combination of them under XOR
	// produces another of the same 16 codewords.
	res.vec = x.vec ^ y.vec
	return res, nil
}

func (cl *CL) setThetaByVec(v1, v2, val uint) error {
	i1, ok := cl.vm[v1]
	if !ok {
		return fmt.Errorf("Vector %x not in vector space", v1)
	}
	i2, ok := cl.vm[v2]
	if !ok {
		return fmt.Errorf("Vector %x not in vector space", v1)
	}
	if i1 > 1<<cl.basisLen || i2 > 1<<cl.basisLen {
		return fmt.Errorf("[!Internal!] Args to setThetaByVec (%x, %x) overflow bitstring of len %d", v1, v2, cl.size*cl.size)
	}
	if val > 0 {
		cl.theta.SetBit(int(i1<<cl.basisLen | i2))
	}
	return nil
}

func (cl *CL) ThetaByVec(v1, v2 uint) (uint, error) {
	if cl.theta == nil {
		return 0, fmt.Errorf("Theta not initialized. Call CL.BuildTheta() first.")
	}

	i1, ok := cl.vm[v1]
	if !ok {
		return 0, fmt.Errorf("Vector %x not in vector space", v1)
	}
	i2, ok := cl.vm[v2]
	if !ok {
		return 0, fmt.Errorf("Vector %x not in vector space", v1)
	}
	if i1 > 1<<cl.basisLen || i2 > 1<<cl.basisLen {
		return 0, fmt.Errorf("[!Internal!] Args to ThetaByVec (%x, %x) overflow bitstring of len %d", v1, v2, cl.size*cl.size)
	}
	return uint(cl.theta.GetBit(int(i1<<cl.basisLen | i2))), nil
}

func (cl *CL) setThetaByIdx(i1, i2, val uint) error {
	if i1 > 1<<cl.basisLen || i2 > 1<<cl.basisLen {
		return fmt.Errorf("Args to setThetaByIdx (%x, %x) overflow bitstring of len %d", i1, i2, cl.size*cl.size)
	}
	if val > 0 {
		cl.theta.SetBit(int(i1<<cl.basisLen | i2))
	}
	return nil
}

func (cl *CL) ThetaByIdx(i1, i2 uint) (uint, error) {

	if cl.theta == nil {
		return 0, fmt.Errorf("Theta not initialized. Call CL.BuildTheta() first.")
	}

	if i1 > 1<<cl.basisLen || i2 > 1<<cl.basisLen {
		return 0, fmt.Errorf("Args to ThetaByIdx (%x, %x) overflow bitstring of len %d", i1, i2, cl.size*cl.size)
	}
	return uint(cl.theta.GetBit(int(i1<<cl.basisLen | i2))), nil
}

func (cl *CL) BuildTheta() error {

	// cf Griess Jr, Robert L. "Code loops." (1986), 226-7

	b0 := cl.basis[0]
	// We build a chain of subspaces V0<V1<V2...<Vk, and define Wk to be
	// Vk+1 - Vk.
	// Intuitively, Wk contains the vectors that will be added when we add the
	// next basis vector bk

	// basic assumption about V0
	e := cl.setThetaByVec(b0, b0, (BitWeight(b0)/4)%2)
	if e != nil {
		return e
	}

	Vk := []uint{0, b0}

	for _, bk := range cl.basis[1:] {

		// create Wk by combining bk with every vector in Vk
		Wk := []uint{}
		for _, v := range Vk {
			Wk = append(Wk, bk^v)
		}

		// D1 - deduce {bk} x Vk and Vk x {bk}
		for _, v := range Vk {
			// Implicitly allow theta(bk, v) to stay set to 0, which is an arbitrary choice.
			cl.setThetaByVec(v, bk, (BitWeight(v&bk)/2)%2)
		}

		// D2 - deduce {bk} x Wk and Wk x {bk}
		for _, v := range Vk {
			a, e := cl.ThetaByVec(bk, v)
			if e != nil {
				return e
			}
			// It looks weird that we're looping over Vk here, but remember
			// that bk^v is an element of _Wk_ not Vk
			cl.setThetaByVec(bk, bk^v, (BitWeight(bk)/4+uint(a))%2)
			cl.setThetaByVec(bk^v, bk, (BitWeight(bk&(bk^v))/2+(BitWeight(bk)/4+uint(a))%2)%2)
		}

		// D3 - deduce Wk x Wk
		for _, v := range Vk {
			for _, w := range Wk {
				a, e := cl.ThetaByVec(v, bk)
				if e != nil {
					return e
				}
				b, e := cl.ThetaByVec(v, bk^w)
				if e != nil {
					return e
				}
				c, e := cl.ThetaByVec(w, bk)
				if e != nil {
					return e
				}
				cl.setThetaByVec(w, bk^v, (BitWeight(v&w)/2+uint(a)+uint(b)+uint(c))%2)
			}
		}

		// D4 - deduce Wk x Vk and Vk x Wk
		for _, w := range Wk {
			for _, v := range Vk {
				a, e := cl.ThetaByVec(w, v^w)
				if e != nil {
					return e
				}
				cl.setThetaByVec(w, v, (BitWeight(w)/4+uint(a))%2)
				cl.setThetaByVec(v, w, (BitWeight(v&w)/2+(BitWeight(w)/4+uint(a))%2)%2)
			}
		}

		Vk = append(Vk, Wk...)

	}

	return nil
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

// Older and slower, but uses full Mul() calls. Useful in case of bugs.
func moufangParallelWorker(work *WorkUnit) {

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
		work.BaseCL.Mul(z, y, zy)
		work.BaseCL.Mul(x, zy, x_zy)
		work.BaseCL.Mul(z, x_zy, z__x_zy)

		//RHS 1
		work.BaseCL.Mul(z, x, zx)
		work.BaseCL.Mul(zx, z, zx_z)
		work.BaseCL.Mul(zx_z, y, zx_z__y)

		// The vector ops are associative, the only thing that can mismatch is the sign.
		fmt.Printf("%s %s %s --> %s vs %s\n", x, y, z, z__x_zy, zx_z__y)

		if z__x_zy.sgn != zx_z__y.sgn {
			work.Failures <- []*CLElem{x, y, z}
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
	elemSz := cl.basisLen + 1
	lim := 1 << (elemSz * 3)
	step := lim / runtime.NumCPU() // more or less
	for cpu := 0; cpu < 1; cpu++ {
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
