package codeloops

import (
	"fmt"
	"github.com/nigelredding/BitString"
	"math/rand"
	"time"
)

// [Gri86] - Griess Jr, Robert L. "Code loops." (1986)

// CL is a structure which contains information about the ambient loop space
// such as the basis.
type CL struct {
	basis     []uint
	basisLen  uint // elems in the basis
	size      uint // elems in the loop
	theta     *Bitstring.Bitstring
	ThetaPath string
	vs        []uint
	vm        map[uint]uint
}

type CLParams struct {
	Basis []uint
	Theta uint
	Seed  int64
}

// NewCL returns a new code loop from a basis for a doubly even binary code.
func NewCL(p CLParams) (cl *CL, e error) {
	cl = new(CL)
	cl.basis = p.Basis
	cl.basisLen = uint(len(p.Basis))
	cl.size = 1 << cl.basisLen
	cl.vs = cl.VectorSpace()
	cl.vm = cl.VectorIdxMap()
	cl.theta = Bitstring.NewBitstring(int(cl.size * cl.size))
	e = cl.buildTheta(p.Theta, p.Seed)
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
		// build the whole vector space by taking linear combinations of all
		// the basis vectors
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

	// This is handy if you want to quickly find out the index for a vector.

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
		return fmt.Errorf("Vector %x not in vector space", v1)
	}
	if i1 >= 1<<cl.basisLen || i2 >= 1<<cl.basisLen {
		return fmt.Errorf("Args to setThetaByVec (%x, %x) overflow bitstring of len %d", v1, v2, cl.size*cl.size)
	}
	if val > 0 {
		cl.theta.SetBit(int(i1<<cl.basisLen | i2))
	}
	return nil
}

func (cl *CL) ThetaByVec(v1, v2 uint) (uint, error) {
	i1, ok := cl.vm[v1]
	if !ok {
		return 0, fmt.Errorf("Vector %x not in vector space", v1)
	}
	i2, ok := cl.vm[v2]
	if !ok {
		return 0, fmt.Errorf("Vector %x not in vector space", v1)
	}
	return uint(cl.theta.GetBit(int(i1<<cl.basisLen | i2))), nil
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

func (cl *CL) ThetaByIdx(i1, i2 uint) (uint, error) {
	if i1 >= 1<<cl.basisLen || i2 >= 1<<cl.basisLen {
		return 0, fmt.Errorf("Args to ThetaByIdx (%x, %x) overflow bitstring of len %d", i1, i2, cl.size*cl.size)
	}
	return uint(cl.theta.GetBit(int(i1<<cl.basisLen | i2))), nil
}

func (cl *CL) thetaByIdxFast(i1, i2 uint) uint {
	return uint(cl.theta.GetBit(int(i1<<cl.basisLen | i2)))
}

func (cl *CL) buildTheta(pathGiven uint, seed int64) error {

	// cf [Gri86] 226-7
	i := uint(0)
	random := false
	if pathGiven == RandomTheta {
		random = true
		if seed == 0 {
			rand.Seed(time.Now().UnixNano())
		} else {
			rand.Seed(seed)
		}
	}
	// we might not take exactly the given path. If the caller supplies a 1
	// for a position that must be 0, we ignore it.
	var pathTaken uint

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
		for _, v := range Vk {
			var x uint
			// theta(bk,0) must be 0 (normalized cocycle), but anything else is up for grabs.
			if v != 0 {
				if random {
					x = uint(rand.Uint32() & 1) // a random bit
				} else {
					x = pathGiven & 1
				}
				cl.setThetaByVec(bk, v, x)
				cl.setThetaByVec(v, bk, ((BitWeight(v&bk)/2)+x)%2)
			} else {
				// theta(bk,v) is implicitly 'set' to 0 in the bitstring
				cl.setThetaByVec(v, bk, (BitWeight(v&bk)/2)%2) // x is forced to be 0
			}
			pathGiven >>= 1
			pathTaken |= x << i
			i++
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
	cl.ThetaPath = fmt.Sprintf("0x%x", pathTaken)
	return nil
}
