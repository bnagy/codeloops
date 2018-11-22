package codeloops

import (
	"fmt"
	"github.com/nigelredding/BitString"
	// "log"
	"math/rand"
	"strconv"
	"strings"
	"time"
)

// [Gri86] - Griess Jr, Robert L. "Code loops." (1986)

// CL is a structure which contains information about the ambient loop space
// such as the basis.
type CL struct {
	basis    []uint
	params   CLParams
	basisLen uint // elems in the basis
	size     uint // elems in the loop
	theta    *Bitstring.Bitstring
	alpha    *Bitstring.Bitstring
	alphaSz  uint
	Seed     string
	vs       []uint
	vsV      []uint
	vsW      []uint
	vsAlpha  []uint
	halfMask uint
	vm       map[uint]uint
	vmAlpha  map[uint]uint
}

type CLParams struct {
	Basis  []uint
	Random bool
	Seed   int64
}

// NewCL returns a new code loop from a basis for a doubly even binary code.
func NewCL(p CLParams) (cl *CL, e error) {
	cl = new(CL)
	cl.params = p
	cl.basis = p.Basis
	cl.basisLen = uint(len(p.Basis))
	cl.size = 1 << cl.basisLen
	cl.vs = VectorSpace(p.Basis)
	m := make(map[uint]uint)
	for i, v := range cl.vs {
		m[v] = uint(i)
	}
	cl.vm = m

	cl.theta = Bitstring.NewBitstring(int(cl.size * cl.size))
	e = cl.buildTheta(p.Random, p.Seed)
	if e != nil {
		return
	}

	// No idea what will happen for odd-length bases
	cl.vsV = VectorSpace(p.Basis[:cl.basisLen/2])
	cl.vsW = VectorSpace(p.Basis[cl.basisLen/2:])
	// HACK this is bigger than it needs to be, because we have a second copy
	// of the zero vector in vsW. The trouble is that all the calculations for
	// setting indicies in the Bitstring work much better if we keep all of
	// the shifts as integers, ie we have the sides of the alpha square stay a
	// power of two.
	cl.vsAlpha = append(cl.vsV, cl.vsW...)
	m2 := make(map[uint]uint)
	for i, v := range cl.vsAlpha {
		m2[v] = uint(i)
	}
	// HACK make sure the index for the zero vector is still zero, since it
	// was written twice above.
	m2[0] = 0
	cl.vmAlpha = m2
	mask, _ := strconv.ParseUint(strings.Repeat("1", int(cl.basisLen/2)), 2, 0)
	cl.halfMask = uint(mask)
	cl.alphaSz = uint(len(cl.vsV) + len(cl.vsW))
	cl.alpha = Bitstring.NewBitstring(int(cl.alphaSz * cl.alphaSz))
	e = cl.buildAlpha()

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
	// code space. This builds any vector as a linear combinations of basis vectors.
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

// Mul performs "multiplication" (the loop action) in the loop.
func (cl *CL) Mul(x, y, res *CLElem) (*CLElem, error) {
	// theta is the cocycle that maps from C^2 -> {0,1}
	t, err := cl.ThetaByVec(x.vec, y.vec)
	if err != nil {
		return nil, err
	}
	// xor is addition in the ambient vector field.
	res.sgn = x.sgn ^ y.sgn ^ t
	res.vec = x.vec ^ y.vec
	return res, nil
}

// MulAlpha performs multiplication, but gets theta values by calculation
// using the much smaller Alpha square.
func (cl *CL) MulAlpha(x, y, res *CLElem) (*CLElem, error) {
	// theta is the cocycle that maps from C^2 -> {0,1}
	t, err := cl.ThetaAlphaByVec(x.vec, y.vec)
	if err != nil {
		return nil, err
	}
	// xor is addition in the ambient vector field.
	res.sgn = x.sgn ^ y.sgn ^ t
	res.vec = x.vec ^ y.vec
	return res, nil
}

// Size returns the number of elements in the loop.
func (cl *CL) Size() int {
	return int(cl.size)
}

// LoopElems returns the loop elements as a slice, with all positive elements
// listed first.
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

// VectorSpace returns a slice of the underlying vectors. For the code loops
// we are working with, these are just the unsigned loop elements.
func (cl *CL) VectorSpace() (vecs []uint) {
	return cl.vs
}

// VectorIdxMap returns a map from vectors to their index in the (ordered) VectorSpace.
func (cl *CL) VectorIdxMap() (m map[uint]uint) {
	// This is handy if you want to quickly find out the index for a vector.
	return cl.vm
}

// Decompose splits a vector into two vectors that can be expressed as
// respective linear combinations of the first and last half of the basis.
func (cl *CL) Decompose(vec uint) (v, w uint, e error) {
	// We build vectors by taking an index 0 <= i < 4096, and then XORing in
	// each basis vector which corresponds to a 1 bit in that index. That
	// means that if we split the index of any vector into the top and bottom
	// 6 bits, we can treat each half as an index into a vector space of size
	// 64 (2^6).
	idx := cl.vm[vec]
	lower := uint(idx) & cl.halfMask
	upper := idx >> (cl.basisLen / 2)
	v = cl.vsV[lower]
	w = cl.vsW[upper]
	if v^w != vec {
		return 0, 0, fmt.Errorf("Failed to decompose vec 0x%x", vec)
	}
	return
}

// PrintBasis displays the loop basis.
func (cl *CL) PrintBasis() {
	fmt.Print("----------\n")
	for i, vec := range cl.basis {
		fmt.Printf("%.2d: %.8b (%d) \n", i, vec, vec)
	}
	fmt.Print("----------\n")
}

// PrintLoopElems displays the loop elements.
func (cl *CL) PrintLoopElems() {
	fmt.Print("----------\n")
	for i, cle := range cl.LoopElems() {
		fmt.Printf("%.2d: %.8b (%s) \n", i, cle.vec, cle.String())
	}
	fmt.Print("----------\n")
}

// PrintVectorSpace displays the underlying vector space.
func (cl *CL) PrintVectorSpace() {
	fmt.Print("----------\n")
	for i, vec := range cl.VectorSpace() {
		fmt.Printf("%.2d: %.8b (%x) \n", i, vec, vec)
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

func (cl *CL) verifyMoufang() (e error) {

	// This version is slow, but works from the basic definition of a Moufang Loop

	elems := cl.LoopElems()
	x, y, z := new(CLElem), new(CLElem), new(CLElem)
	zy, x_zy, z__x_zy := new(CLElem), new(CLElem), new(CLElem)
	zx, zx_z, zx_z__y := new(CLElem), new(CLElem), new(CLElem)

	for i := 0; i < len(elems); i++ {
		for j := 0; j < len(elems); j++ {
			for k := 0; k < len(elems); k++ {

				x, y, z = &elems[i], &elems[j], &elems[k]

				// LHS
				cl.Mul(z, y, zy)
				cl.Mul(x, zy, x_zy)
				cl.Mul(z, x_zy, z__x_zy)

				//RHS
				cl.Mul(z, x, zx)
				cl.Mul(zx, z, zx_z)
				cl.Mul(zx_z, y, zx_z__y)

				// The vector ops are associative, the only thing that can mismatch is the sign.
				if z__x_zy.sgn != zx_z__y.sgn {
					return fmt.Errorf("Code loop failed Moufang identity at %s, %s, %s", x, y, z)
				}
			}
		}
	}

	return nil
}

func (cl *CL) verifyMoufang2() (e error) {

	// This version is an order of magnitude faster, and uses an identity from [Gri86] p. 226
	// Moufang iff t(x,y) + t (z,x) + t(x+y, z+x) = t(y,z) + t(x, y+z) + t(x+y+z, x)

	vs := cl.VectorSpace()
	var x, y, z uint
	for i := 0; i < len(vs); i++ {
		for j := 0; j < len(vs); j++ {
			for k := 0; k < len(vs); k++ {
				x, y, z = vs[i], vs[j], vs[k]
				if (cl.thetaByVecFast(x, y)^
					cl.thetaByVecFast(z, x)^
					cl.thetaByVecFast(x^y, z^x))^
					cl.thetaByVecFast(y, z)^
					cl.thetaByVecFast(x, y^z)^
					cl.thetaByVecFast(x^y^z, x) != 0 {
					return fmt.Errorf("Code loop failed Moufang identity at %x, %x, %x", x, y, z)

				}
			}
		}
	}
	return nil
}

// IsMoufang checks whether the loop is Moufang.
func (cl *CL) IsMoufang() bool {
	e := cl.verifyMoufang2()
	return e == nil
}

func (cl *CL) verifyAssoc() (isAssoc bool) {

	// Traditional check - for all a,b,c (ab)c = a(bc)
	// (this runs on the loop elems)

	isAssoc = true
	elems := cl.LoopElems()

	SetCombinationsWithReplacement(uint(len(elems)), 3, func(triple []uint) error {
		a, b, c := &elems[triple[0]], &elems[triple[1]], &elems[triple[2]]
		a_bc, ab_c := new(CLElem), new(CLElem)

		cl.Mul(b, c, a_bc)
		cl.Mul(a, a_bc, a_bc)
		cl.Mul(a, b, ab_c)
		cl.Mul(ab_c, c, ab_c)

		if a_bc.sgn != ab_c.sgn {
			isAssoc = false
			return fmt.Errorf("Not assoc")
		}
		return nil
	})
	return
}

func (cl *CL) verifyAssoc2() bool {

	// Faster check, using the cocycle identity
	// (this only runs on the vector space elems; half the size)

	vs := cl.VectorSpace()
	var x, y, z uint
	for i := 0; i < len(vs); i++ {
		for j := 0; j < len(vs); j++ {
			for k := 0; k < len(vs); k++ {
				x, y, z = vs[i], vs[j], vs[k]
				if (cl.thetaByVecFast(x, y^z) ^ cl.thetaByVecFast(y, z) ^
					cl.thetaByVecFast(x^y, z) ^ cl.thetaByVecFast(x, y)) != 0 {
					return false

				}
			}
		}
	}
	return true
}

// IsAssoc checks whether the loop is fully associative (ie a group).
func (cl *CL) IsAssoc() (isAssoc bool) {
	return cl.verifyAssoc2()
}

func (cl *CL) setThetaByVec(v1, v2, val uint) error {
	i1, ok := cl.vm[v1]
	if !ok {
		return fmt.Errorf("Vector %x not in vector space", v1)
	}
	i2, ok := cl.vm[v2]
	if !ok {
		return fmt.Errorf("Vector %x not in vector space", v2)
	}
	if i1 >= 1<<cl.basisLen || i2 >= 1<<cl.basisLen {
		return fmt.Errorf("Args to setThetaByVec (%x, %x) overflow bitstring of len %d", v1, v2, cl.size*cl.size)
	}
	if val > 0 {
		cl.theta.SetBit(int(i1<<cl.basisLen | i2))
	}
	return nil
}

// ThetaByVec returns theta(v1, v2) where v1 and v2 are vectors in the underlying space.
func (cl *CL) ThetaByVec(v1, v2 uint) (uint, error) {
	i1, ok := cl.vm[v1]
	if !ok {
		return 0, fmt.Errorf("Vector %x not in vector space", v1)
	}
	i2, ok := cl.vm[v2]
	if !ok {
		return 0, fmt.Errorf("Vector %x not in vector space", v2)
	}
	return uint(cl.theta.GetBit(int(i1<<cl.basisLen | i2))), nil
}

func (cl *CL) ThetaAlphaByVec(x, y uint) (uint, error) {
	v1, w1, err := cl.Decompose(x)
	if err != nil {
		return 0, err
	}
	v2, w2, err := cl.Decompose(y)
	if err != nil {
		return 0, err
	}
	// With the split basis this will always be 0
	a, err := cl.AlphaByVec(v1, v2)
	if err != nil {
		return 0, err
	}
	b, err := cl.AlphaByVec(w1, w2)
	if err != nil {
		return 0, err
	}
	c, err := cl.AlphaByVec(v1, w1)
	if err != nil {
		return 0, err
	}
	d, err := cl.AlphaByVec(w2, v2)
	if err != nil {
		return 0, err
	}
	e, err := cl.AlphaByVec(v1^v2, w1^w2)
	if err != nil {
		return 0, err
	}
	f := BitWeight(v2&(w1^w2)) / 2
	g := BitWeight(v1 & v2 & (w1 ^ w2))
	h := BitWeight((w1 & w2 & v2))
	i := BitWeight(v1 & w1 & (v2 ^ w2))
	return (a + b + c + d + e + f + g + h + i) % 2, nil
}

// Despite the name, this is not actually much faster.
func (cl *CL) thetaAlphaByVecFast(x, y uint) uint {
	v1, w1, _ := cl.Decompose(x)
	v2, w2, _ := cl.Decompose(y)
	a := cl.alphaByVecFast(v1, v2)
	b := cl.alphaByVecFast(w1, w2)
	c := cl.alphaByVecFast(v1, w1)
	d := cl.alphaByVecFast(w2, v2)
	e := cl.alphaByVecFast(v1^v2, w1^w2)
	f := BitWeight(v2&(w1^w2)) / 2
	g := BitWeight(v1 & v2 & (w1 ^ w2))
	h := BitWeight((w1 & w2 & v2))
	i := BitWeight(v1 & w1 & (v2 ^ w2))
	return (a + b + c + d + e + f + g + h + i) % 2
}

func (cl *CL) thetaByVecFast(v1, v2 uint) uint {
	return uint(cl.theta.GetBit(int(cl.vm[v1]<<cl.basisLen | cl.vm[v2])))
}

func (cl *CL) setThetaByIdx(i1, i2, val uint) error {
	if i1 >= 1<<cl.basisLen || i2 >= 1<<cl.basisLen {
		return fmt.Errorf("Args to setThetaByIdx (%x, %x) overflow bitstring of len %d", i1, i2, cl.size*cl.size)
	}
	if val > 0 {
		cl.theta.SetBit(int(i1<<cl.basisLen | i2))
	}
	return nil
}

// ThetaByIdx returns theta(i1, i2) where i1 and i2 are indicies into the vector space.
func (cl *CL) ThetaByIdx(i1, i2 uint) (uint, error) {
	if i1 >= 1<<cl.basisLen || i2 >= 1<<cl.basisLen {
		return 0, fmt.Errorf("Args to ThetaByIdx (%x, %x) overflow bitstring of len %d", i1, i2, cl.size*cl.size)
	}
	return uint(cl.theta.GetBit(int(i1<<cl.basisLen | i2))), nil
}

func (cl *CL) thetaByIdxFast(i1, i2 uint) uint {
	return uint(cl.theta.GetBit(int(i1<<cl.basisLen | i2)))
}

func (cl *CL) buildTheta(random bool, seed int64) error {

	var x uint

	if random && seed == 0 {
		rand.Seed(time.Now().UnixNano())
	} else {
		rand.Seed(seed)
	}

	// cf [Gri86] 226-7
	// We build a chain of subspaces V0<V1<V2...<Vk, and define Wk to be
	// Vk+1 - Vk.
	// Intuitively, Wk contains the vectors that will be added when we add the
	// next basis vector bk

	b0 := cl.basis[0]
	// basic assumption about V0, other than normalization, which is 'set' for
	// free, since the bitstring defaults to 00..0
	e := cl.setThetaByVec(b0, b0, (BitWeight(b0)/4)%2)
	if e != nil {
		return e
	}

	Vk := []uint{0, b0} // span(b0)

	for _, bk := range cl.basis[1:] {

		// create Wk by combining bk with every vector in Vk
		Wk := []uint{}
		for _, v := range Vk {
			Wk = append(Wk, bk^v)
		}

		// D1 - define {bk} x Vk and deduce Vk x {bk}
		// log.Printf("Starting D1")
		for _, v := range Vk {
			// theta(bk,0) must be 0 (normalized cocycle), but anything else is up for grabs.
			if v != 0 {
				if random {
					x = uint(rand.Uint32() & 1) // a random bit
				} else {
					x = 0
				}
				cl.setThetaByVec(bk, v, x)
				cl.setThetaByVec(v, bk, ((BitWeight(v&bk)/2)+x)%2)
			} else {
				// theta(bk,v) is implicitly 'set' to 0 in the bitstring
				cl.setThetaByVec(v, bk, (BitWeight(v&bk)/2)%2) // x is forced to be 0
			}
		}

		// D2 - deduce {bk} x Wk and Wk x {bk}
		// log.Printf("Starting D2")
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
		// log.Printf("Starting D3")
		for _, v := range Vk {
			for _, v2 := range Vk {
				w := bk ^ v2
				a, e := cl.ThetaByVec(v, bk)
				if e != nil {
					e := fmt.Errorf("Error getting %.2x, %.2x: %s", v, bk)
					return e
				}
				b, e := cl.ThetaByVec(v, bk^w)
				if e != nil {
					e := fmt.Errorf("Error getting %.2x, %.2x: %s", v, bk^w)
					return e
				}
				c, e := cl.ThetaByVec(w, bk)
				if e != nil {
					e := fmt.Errorf("Error getting %.2x, %.2x: %s", w, bk)
					return e
				}
				res := (BitWeight(v&w)/2 + uint(a) + uint(b) + uint(c)) % 2
				cl.setThetaByVec(w, bk^v, res)
			}
		}

		// D4 - deduce Wk x Vk and Vk x Wk
		// log.Printf("Starting D4")
		for _, v := range Vk {
			for _, v2 := range Vk {
				w := bk ^ v2
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
	cl.Seed = fmt.Sprintf("0x%x", seed)
	return nil
}

func (cl *CL) setAlphaByVec(v1, v2, val uint) error {

	i1, ok := cl.vmAlpha[v1]
	if !ok {
		return fmt.Errorf("v1 %x not in vector space", v1)
	}
	i2, ok := cl.vmAlpha[v2]
	if !ok {
		return fmt.Errorf("v2 %x not in vector space", v2)
	}
	if i1 >= cl.alphaSz || i2 >= cl.alphaSz {
		return fmt.Errorf("Args to setAlphaByVec (%x, %x) overflow bitstring of len %d", v1, v2, cl.alphaSz*cl.alphaSz)
	}

	// to optimise, we could check this right at the start, but I want the
	// error checking to run.
	if val > 0 {
		cl.alpha.SetBit(int((i1 * cl.alphaSz) | i2))
	}
	return nil
}

func (cl *CL) buildAlpha() (e error) {
	for i := uint(0); i < cl.alphaSz; i++ {
		for j := uint(0); j < cl.alphaSz; j++ {
			res, err := cl.ThetaByVec(cl.vsAlpha[i], cl.vsAlpha[j])
			if err != nil {
				e = fmt.Errorf("buildAlpha: %s", e)
				return
			}
			err = cl.setAlphaByVec(cl.vsAlpha[i], cl.vsAlpha[j], res)
			if err != nil {
				e = fmt.Errorf("buildAlpha: %s", err)
				return
			}
		}
	}
	return
}

// AlphaByVec returns alpha(v1, v2) where v1 and v2 are vectors in the Alpha square
func (cl *CL) AlphaByVec(v1, v2 uint) (uint, error) {
	i1, ok := cl.vmAlpha[v1]
	if !ok {
		return 0, fmt.Errorf("Vector %x not in alpha space", v1)
	}
	i2, ok := cl.vmAlpha[v2]
	if !ok {
		return 0, fmt.Errorf("Vector %x not in alpha space", v2)
	}
	return uint(cl.alpha.GetBit(int(i1*cl.alphaSz | i2))), nil
}

func (cl *CL) alphaByVecFast(v1, v2 uint) uint {
	return uint(cl.alpha.GetBit(int(cl.vmAlpha[v1]*cl.alphaSz | cl.vmAlpha[v2])))
}

// AlphaByIdx returns alpha(i1, i2) where i1 and i2 are indicies into the Alpha square.
func (cl *CL) AlphaByIdx(i1, i2 uint) (uint, error) {
	if i1 >= cl.alphaSz || i2 >= cl.alphaSz {
		return 0, fmt.Errorf("Args to AlphaByIdx (%x, %x) overflow bitstring of len %d", i1, i2, cl.alphaSz*cl.alphaSz)
	}
	return uint(cl.alpha.GetBit(int(i1*cl.alphaSz | i2))), nil
}

// This is legacy code which is here in case I ever need to regenerate some
// old images. The issue is that the path supplied is too small, so what's
// actually happening is that only the first 2^64 choices are being changed by
// a given path and the rest ends up being deterministic.
// func (cl *CL) buildThetaWithPath(pathGiven uint, seed int64) error {

// 	var x uint

// 	// cf [Gri86] 226-7
// 	i := uint(0)
// 	random := false
// 	if pathGiven == RandomTheta {
// 		random = true
// 		if seed == 0 {
// 			rand.Seed(time.Now().UnixNano())
// 		} else {
// 			rand.Seed(seed)
// 		}
// 	}
// 	// we might not take exactly the given path. If the caller supplies a 1
// 	// for a position that must be 0, we ignore it.
// 	var pathTaken uint

// 	// We build a chain of subspaces V0<V1<V2...<Vk, and define Wk to be
// 	// Vk+1 - Vk.
// 	// Intuitively, Wk contains the vectors that will be added when we add the
// 	// next basis vector bk

// 	b0 := cl.basis[0]
// 	// basic assumption about V0, other than normalization, which is 'set' for
// 	// free, since the bitstring defaults to 00..0
// 	e := cl.setThetaByVec(b0, b0, (BitWeight(b0)/4)%2)
// 	if e != nil {
// 		return e
// 	}

// 	Vk := []uint{0, b0} // span(b0)

// 	for _, bk := range cl.basis[1:] {

// 		// create Wk by combining bk with every vector in Vk
// 		Wk := []uint{}
// 		for _, v := range Vk {
// 			Wk = append(Wk, bk^v)
// 		}

// 		// D1 - define {bk} x Vk and deduce Vk x {bk}
// 		// log.Printf("Starting D1")
// 		for _, v := range Vk {
// 			// theta(bk,0) must be 0 (normalized cocycle), but anything else is up for grabs.
// 			if v != 0 {
// 				if random {
// 					x = uint(rand.Uint32() & 1) // a random bit
// 				} else {
// 					x = pathGiven & 1
// 				}
// 				cl.setThetaByVec(bk, v, x)
// 				cl.setThetaByVec(v, bk, ((BitWeight(v&bk)/2)+x)%2)
// 			} else {
// 				// theta(bk,v) is implicitly 'set' to 0 in the bitstring
// 				cl.setThetaByVec(v, bk, (BitWeight(v&bk)/2)%2) // x is forced to be 0
// 			}
// 			pathGiven >>= 1
// 			pathTaken |= x << i
// 			i++
// 		}

// 		// D2 - deduce {bk} x Wk and Wk x {bk}
// 		// log.Printf("Starting D2")
// 		for _, v := range Vk {
// 			a, e := cl.ThetaByVec(bk, v)
// 			if e != nil {
// 				return e
// 			}
// 			// It looks weird that we're looping over Vk here, but remember
// 			// that bk^v is an element of _Wk_ not Vk
// 			cl.setThetaByVec(bk, bk^v, (BitWeight(bk)/4+uint(a))%2)
// 			cl.setThetaByVec(bk^v, bk, (BitWeight(bk&(bk^v))/2+(BitWeight(bk)/4+uint(a))%2)%2)
// 		}

// 		// D3 - deduce Wk x Wk
// 		// log.Printf("Starting D3")
// 		for _, v := range Vk {
// 			for _, v2 := range Vk {
// 				w := bk ^ v2
// 				a, e := cl.ThetaByVec(v, bk)
// 				if e != nil {
// 					e := fmt.Errorf("Error getting %.2x, %.2x: %s", v, bk)
// 					return e
// 				}
// 				b, e := cl.ThetaByVec(v, bk^w)
// 				if e != nil {
// 					e := fmt.Errorf("Error getting %.2x, %.2x: %s", v, bk^w)
// 					return e
// 				}
// 				c, e := cl.ThetaByVec(w, bk)
// 				if e != nil {
// 					e := fmt.Errorf("Error getting %.2x, %.2x: %s", w, bk)
// 					return e
// 				}
// 				res := (BitWeight(v&w)/2 + uint(a) + uint(b) + uint(c)) % 2
// 				cl.setThetaByVec(w, bk^v, res)
// 			}
// 		}

// 		// D4 - deduce Wk x Vk and Vk x Wk
// 		// log.Printf("Starting D4")
// 		for _, v := range Vk {
// 			for _, v2 := range Vk {
// 				w := bk ^ v2
// 				a, e := cl.ThetaByVec(w, v^w)
// 				if e != nil {
// 					return e
// 				}
// 				cl.setThetaByVec(w, v, (BitWeight(w)/4+uint(a))%2)
// 				cl.setThetaByVec(v, w, (BitWeight(v&w)/2+(BitWeight(w)/4+uint(a))%2)%2)
// 			}
// 		}

// 		Vk = append(Vk, Wk...)

// 	}
// 	cl.Seed = fmt.Sprintf("0x%x", pathTaken)
// 	return nil
// }
