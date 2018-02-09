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

func TestVerifyHammingMoufangGood(t *testing.T) {
	cl, err := NewCL(HammingBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.VerifyMoufang()
	if err != nil {
		t.Fatalf("Hamming basis failed VerifyMoufang(): %s", err)
	}
}

func TestVerifyHammingNotAssoc(t *testing.T) {
	cl, err := NewCL(GolayBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	isAssoc := true
	elems := cl.LoopElems()
	f := func(triple []uint) error {
		a, b, c := &elems[triple[0]], &elems[triple[1]], &elems[triple[2]]
		// fmt.Printf("a: %s b: %s c: %s ", a, b, c)

		a_bc, ab_c := new(CLElem), new(CLElem)
		a_bc.Mul(b, c)
		// fmt.Printf("bc: %s ", a_bc)
		a_bc.Mul(a, a_bc)
		ab_c.Mul(a, b)
		// fmt.Printf("ab: %s ", ab_c)
		ab_c.Mul(ab_c, c)
		// fmt.Printf("a(bc): %s (ab)c: %s\n", a_bc, ab_c)

		if a_bc.sgn != ab_c.sgn {
			isAssoc = false
			return fmt.Errorf("Not assoc")
		}
		return nil
	}
	SetCombinationsWithReplacement(uint(len(elems)), 3, f)
	if isAssoc {
		t.Fatalf("Hamming basis produced an associative group, expected a loop.")
	}
}

func TestVerifyHammingMoufangBad(t *testing.T) {
	cl, err := NewCL(badHammingBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.VerifyMoufang()
	if err == nil {
		t.Fatalf("Bad Hamming basis passed VerifyMoufang(), expected it to fail.")
	}
}

// func TestVerifyGolayMoufangGood(t *testing.T) {
// 	cl, err := NewCL(GolayBasis)
// 	if err != nil {
// 		t.Fatalf("Failed to create CL: %s", err)
// 	}
// 	err = cl.VerifyMoufang()
// 	if err != nil {
// 		t.Fatalf("Golay basis failed VerifyMoufang(): %s", err)
// 	}
// }

func TestVerifyGolayMoufangBad(t *testing.T) {
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

func TestBuildTheta(t *testing.T) {
	basis := HammingBasis
	cl, err := NewCL(basis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	theta, err := cl.buildTheta2()
	if err != nil {
		t.Fatal(err)
	}
	if len(theta) != 1<<(uint(len(basis)*2)) {
		fmt.Printf("%#v\n", cl.VectorSpace())

		for k, v := range theta {
			fmt.Printf("%.4x %d\n", k, v)
		}

		t.Fatalf("Theta map incomplete. Expected %d entries, got %d", 1<<(uint(len(basis)*2)), len(theta))
	}

	vs := cl.VectorSpace()

	for k, v := range theta {
		fmt.Printf("%.4x %d\n", k, v)
	}
	// Test theta properties.

	// theta(0,x) == theta(x,0) == 0 (normalized cocycle)
	for k, v := range theta {
		if k&0xff == 0 || k&0xff00 == 0 {
			// upper or lower 8 bits are 0, so this entry should be normalized
			if v != 0 {
				t.Fatalf("Theta not normalised - expected 0 for %x, got %d", k, v)
			}

		}
	}

	// for all x, theta(x,x) === |x|/4
	for i := 0; i < len(vs); i++ {
		x := vs[i]
		if _, ok := theta[x<<8|x]; !ok {
			t.Fatalf("Missing entry for theta(%x,%x)", x, x)
		}
		if theta[x<<8|x] != (BitWeight(x)/4)%2 {
			t.Fatalf("Expected theta(%x,%x) to be %d, got %d", x, x, (BitWeight(x)/4)%2, theta[x<<8|x])
		}
	}

	// for all x,y, theta(x,y) + theta(y,x) === |x&y|/2
	for i := 0; i < len(vs); i++ {
		for j := 0; j < len(vs); j++ {
			x, y := vs[i], vs[j]
			if _, ok := theta[x<<8|y]; !ok {
				t.Fatalf("Missing entry for theta(%x,%x)", x, y)
			}
			if _, ok := theta[y<<8|x]; !ok {
				t.Fatalf("Missing entry for theta(%x,%x)", y, x)
			}
			lhs := (theta[x<<8|y] + theta[y<<8|x]) % 2
			rhs := (BitWeight(x&y) / 2) % 2
			if rhs != lhs {
				t.Fatalf("Expected theta(%x,%x) + theta(%x,%x) to be %d, got %d", x, y, y, x, rhs, lhs)
			}
		}
	}

	// for all x,y,z theta(x,y) + theta(x^y,z) + theta(y,z) + theta(x,y^z) === |x&y&z|
	for i := 0; i < len(vs); i++ {
		for j := 0; j < len(vs); j++ {
			for k := 0; k < len(vs); k++ {
				x, y, z := vs[i], vs[j], vs[k]
				if _, ok := theta[x<<8|y]; !ok {
					t.Fatalf("Missing entry for theta(%x,%x)", x, y)
				}
				if _, ok := theta[(x^y)<<8|z]; !ok {
					t.Fatalf("Missing entry for theta(%x,%x)", (x ^ y), z)
				}
				if _, ok := theta[y<<8|z]; !ok {
					t.Fatalf("Missing entry for theta(%x,%x)", y, z)
				}
				if _, ok := theta[x<<8|(y^z)]; !ok {
					t.Fatalf("Missing entry for theta(%x,%x)", x, (y ^ z))
				}
				lhs := (theta[x<<8|y] + theta[(x^y)<<8|z] + theta[y<<8|z] + theta[x<<8|(y^z)]) % 2
				rhs := (BitWeight(x & y & z)) % 2
				if rhs != lhs {
					t.Errorf("(x,y - %x,%x): %d (x^y,z - %x,%x): %d (x,y^z - %x,%x): %d (y,z - %x,%x): %d\n",
						x, y, theta[x<<8|y],
						x^y, z, theta[(x^y)<<8|z],
						x, y^z, theta[y<<8|z],
						y, z, theta[x<<8|(y^z)])
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

func BenchmarkSigma(b *testing.B) {
	dummy := uint(0)
	for i := 0; i < b.N; i++ {
		dummy |= sigma(0xdeadbeef, 0xfacef00d)
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
	x, _ := cl.NewCLElem(5, Pos)
	y, _ := cl.NewCLElem(7, Neg)
	res := new(CLElem)
	for n := 0; n < b.N; n++ {
		res.Mul(x, y)
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
