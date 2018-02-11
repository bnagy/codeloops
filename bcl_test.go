package codeloops

import (
	"fmt"
	"testing"
)

var badHammingBasis = []uint{
	0x87, // 1000 | 0111
	0x4b, // 0100 | 1011
	0x3d, // 0011 | 1101 <-- should be 0x2d
	0x1e, // 0001 | 1110
}

var badGolayBasis = []uint{
	0x8009f1,
	0x4004fa,
	0x20027d,
	0x11093e, // <-- should be 0x10093e
	0x80c9d,
	0x40e4e,
	0x20f25,
	0x10f92,
	0x87c9,
	0x43e6,
	0x2557,
	0x1aab,
}

func TestBitWeight(t *testing.T) {
	tests := []uint{0x0, 0x1, 0x87, 0xff, 0xffffffff, 0xdeadbeef, 0xffffffffdeadbeef}
	expect := []uint{0, 1, 4, 8, 32, 24, 56}
	for idx, x := range tests {
		if BitWeight(x) != expect[idx] {
			t.Errorf("TestBitWeight: Weight of 0x%x should be %d, got %d", x, expect[idx], BitWeight(x))
		}
	}
}

func TestPopCntq(t *testing.T) {
	tests := []uint{0x0, 0x1, 0x87, 0xff, 0xffffffff, 0xdeadbeef, 0xffffffffdeadbeef}
	expect := []uint{0, 1, 4, 8, 32, 24, 56}
	for idx, x := range tests {
		if popCntq(x) != expect[idx] {
			t.Errorf("TestPopCntq: Weight of 0x%x should be %d, got %d", x, expect[idx], popCntq(x))
		}
	}
}

func TestVerifyBasis(t *testing.T) {
	cl, err := NewCL(HammingBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.VerifyBasis()
	if err != nil {
		t.Fatalf("Hamming basis failed VerifyBasis(), expected it to pass.")
	}
	cl, err = NewCL(badHammingBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.VerifyBasis()
	if err == nil {
		t.Fatalf("Bad Hamming basis passed VerifyBasis(), expected it to fail.")
	}
	cl, err = NewCL(GolayBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.VerifyBasis()
	if err != nil {
		t.Fatalf("Golay basis failed VerifyBasis(), expected it to pass.")
	}
	cl, err = NewCL(badGolayBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.VerifyBasis()
	if err == nil {
		t.Fatalf("Bad Golay basis passed VerifyBasis(), expected it to fail.")
	}
}

func TestHammingMoufang(t *testing.T) {
	cl, err := NewCL(HammingBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.VerifyMoufang()
	if err != nil {
		t.Fatalf("Hamming loop not Moufang: %s", err)
	}
}

func TestHammingNotAssoc(t *testing.T) {
	cl, err := NewCL(HammingBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	isAssoc := true
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

	if isAssoc {
		t.Fatalf("Hamming basis produced an associative group, expected a loop.")
	}
}

func TestBadHammingNotMoufang(t *testing.T) {
	cl, err := NewCL(badHammingBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.VerifyMoufang()
	if err == nil {
		t.Fatalf("Bad Hamming basis passed Moufang test.")
	}
}

func TestGolayMoufang(t *testing.T) {
	cl, err := NewCL(GolayBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.VerifyMoufang()
	if err != nil {
		t.Fatalf("Golay basis failed VerifyMoufang(): %s", err)
	}
}

func TestBadGolayNotMoufang(t *testing.T) {
	cl, err := NewCL(badGolayBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.VerifyMoufang()
	if err == nil {
		t.Fatalf("Bad Golay basis passed VerifyMoufang(), expected it to fail.")
	}
}

func TestListHammingElems(t *testing.T) {
	cl, err := NewCL(HammingBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	cl.PrintLoopElems()
}

func TestListHammingVectorSpace(t *testing.T) {
	cl, err := NewCL(HammingBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	cl.PrintVectorSpace()
}

func TestHammingTheta(t *testing.T) {

	basis := HammingBasis
	cl, err := NewCL(basis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}

	// Print out the cocycle, for laughs.
	vs := cl.VectorSpace()
	fmt.Printf(" * |")
	for _, v := range vs {
		fmt.Printf("%.2x|", v)
	}
	fmt.Printf("\n")
	for _, v1 := range vs {
		fmt.Printf("|%.2x|", v1)
		for _, v2 := range vs {
			res, e := cl.ThetaByVec(v1, v2)
			if e != nil {
				t.Fatal(e)
			}
			if res == 1 {
				fmt.Printf("XX|")
			} else {
				fmt.Printf("  |")
			}
		}
		fmt.Printf("\n")
	}

	// Verify that the cocycle is normalized, then check the axioms (S), (C),
	// (A) from [Griess86] p. 225

	// theta(0,x) == theta(x,0) == 0 (normalized cocycle)
	for _, v := range vs {
		if b, _ := cl.ThetaByVec(v, 0); b != 0 {
			t.Fatalf("Theta not normalized at %x, 0", v)
		}
		if b, _ := cl.ThetaByVec(0, v); b != 0 {
			t.Fatalf("Theta not normalized at 0, %x", v)
		}
	}

	// (S) for all x, theta(x,x) === |x|/4 (all congruence is mod 2)
	for i := 0; i < len(vs); i++ {
		if b, _ := cl.ThetaByIdx(uint(i), uint(i)); b != (BitWeight(vs[i])/4)%2 {
			t.Fatalf("Expected theta(%x,%x) to be %d, got %d", i, i, (BitWeight(vs[i])/4)%2, b)
		}
	}

	// (C) for all x,y, theta(x,y) + theta(y,x) === |x&y|/2
	for i := 0; i < len(vs); i++ {
		for j := 0; j < len(vs); j++ {
			x, y := vs[i], vs[j]
			a, _ := cl.ThetaByIdx(uint(i), uint(j))
			b, _ := cl.ThetaByIdx(uint(j), uint(i))
			lhs := uint((a + b) % 2)
			rhs := (BitWeight(x&y) / 2) % 2
			if rhs != lhs {
				t.Fatalf("Expected theta(%x,%x) + theta(%x,%x) to be %d, got %d", x, y, y, x, rhs, lhs)
			}
		}
	}

	// (A) for all x,y,z theta(x,y) + theta(x^y,z) + theta(y,z) + theta(x,y^z) === |x&y&z|
	for i := 0; i < len(vs); i++ {
		for j := 0; j < len(vs); j++ {
			for k := 0; k < len(vs); k++ {
				x, y, z := vs[i], vs[j], vs[k]
				a, _ := cl.ThetaByVec(x, y)
				b, _ := cl.ThetaByVec(x^y, z)
				c, _ := cl.ThetaByVec(y, z)
				d, _ := cl.ThetaByVec(x, y^z)
				lhs := (a + b + c + d) % 2
				rhs := (BitWeight(x & y & z)) % 2
				if rhs != lhs {
					t.Errorf("(x,y - %x,%x): %d (x^y,z - %x,%x): %d (x,y^z - %x,%x): %d (y,z - %x,%x): %d\n",
						x, y, a,
						x^y, z, b,
						x, y^z, c,
						y, z, d)
					t.Fatalf("Error in triple identity for %x %x %x, Expected %d, got %d", x, y, z, rhs, lhs)
				}
			}
		}
	}
}

func BenchmarkBitWeight(b *testing.B) {
	// dummy var stops the compiler from eliminating dead code in the loop
	dummy := uint(0)
	for i := 0; i < b.N; i++ {
		dummy |= BitWeight(0xdeadbeef)
	}
}

// This benches much slower, despite being implemented as a single asm
// instruction? Maybe the trampoline process costs a nanosecond or so...
func BenchmarkPopCntq(b *testing.B) {
	dummy := uint(0)
	for i := 0; i < b.N; i++ {
		dummy |= popCntq(0xdeadbeef)
	}
}

func BenchmarkVerifyBasisHamming(b *testing.B) {
	cl, err := NewCL(HammingBasis)
	if err != nil {
		b.Fatalf("Failed to create CL: %s", err)
	}
	for n := 0; n < b.N; n++ {
		err = cl.VerifyBasis()
		if err != nil {
			b.Fatalf("Hamming basis failed VerifyBasis(), expected it to pass.")
		}
	}
}

func BenchmarkVerifyBasisGolay(b *testing.B) {
	cl, err := NewCL(GolayBasis)
	if err != nil {
		b.Fatalf("Failed to create CL: %s", err)
	}
	for n := 0; n < b.N; n++ {
		err = cl.VerifyBasis()
		if err != nil {
			b.Fatalf("Golay basis failed VerifyBasis(), expected it to pass.")
		}
	}
}

func BenchmarkHammingMoufang(b *testing.B) {
	cl, err := NewCL(HammingBasis)
	if err != nil {
		b.Fatalf("Failed to create CL: %s", err)
	}
	for n := 0; n < b.N; n++ {
		err = cl.VerifyMoufang()
		if err != nil {
			b.Fatalf("%s", err)
		}
	}
}

func BenchmarkMul(b *testing.B) {
	cl, err := NewCL(HammingBasis)
	if err != nil {
		b.Fatalf("Failed to create CL: %s", err)
	}
	x, _ := cl.NewElemFromIdx(5, Pos)
	y, _ := cl.NewElemFromIdx(7, Neg)
	res := new(CLElem)
	for n := 0; n < b.N; n++ {
		cl.Mul(x, y, res)
	}

}

func BenchmarkLoopElemsHamming(b *testing.B) {
	cl, err := NewCL(HammingBasis)
	if err != nil {
		b.Fatalf("Failed to create CL: %s", err)
	}
	for n := 0; n < b.N; n++ {
		cl.LoopElems()
	}
}

func BenchmarkLoopElemsGolay(b *testing.B) {
	cl, err := NewCL(GolayBasis)
	if err != nil {
		b.Fatalf("Failed to create CL: %s", err)
	}
	for n := 0; n < b.N; n++ {
		cl.LoopElems()
	}
}
