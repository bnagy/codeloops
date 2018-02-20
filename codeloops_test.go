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
	cl, err := NewCL(CLParams{Basis: HammingBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.VerifyBasis()
	if err != nil {
		t.Fatalf("Hamming basis failed VerifyBasis(), expected it to pass.")
	}
	cl, err = NewCL(CLParams{Basis: badHammingBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.VerifyBasis()
	if err == nil {
		t.Fatalf("Bad Hamming basis passed VerifyBasis(), expected it to fail.")
	}
	cl, err = NewCL(CLParams{Basis: GolayBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.VerifyBasis()
	if err != nil {
		t.Fatalf("Golay basis failed VerifyBasis(), expected it to pass.")
	}
	cl, err = NewCL(CLParams{Basis: badGolayBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.VerifyBasis()
	if err == nil {
		t.Fatalf("Bad Golay basis passed VerifyBasis(), expected it to fail.")
	}
}

func TestHammingMoufang(t *testing.T) {
	for i := 0; i < 1000; i++ {
		cl, err := NewCL(CLParams{Basis: HammingBasis, Theta: RandomTheta})
		if err != nil {
			t.Fatalf("Failed to create CL: %s", err)
		}
		err = cl.verifyMoufang()
		if err != nil {
			t.Fatalf("Hamming loop not Moufang: %s", err)
		}
	}
}

func TestHammingMoufang2(t *testing.T) {
	for i := 0; i < 1000; i++ {
		cl, err := NewCL(CLParams{Basis: HammingBasis, Theta: RandomTheta})
		if err != nil {
			t.Fatalf("Failed to create CL: %s", err)
		}
		err = cl.verifyMoufang2()
		if err != nil {
			t.Fatalf("Hamming loop not Moufang: %s", err)
		}
	}
}

func TestHammingNotAssoc(t *testing.T) {
	cl, err := NewCL(CLParams{Basis: HammingBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	if cl.IsAssoc() {
		t.Fatalf("Hamming basis produced an associative group, expected a loop.")
	}
}

func TestBadHammingNotMoufang(t *testing.T) {
	cl, err := NewCL(CLParams{Basis: badHammingBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.verifyMoufang()
	if err == nil {
		t.Fatalf("Bad Hamming basis passed Moufang test.")
	}
}

func TestBadHammingNotMoufang2(t *testing.T) {
	cl, err := NewCL(CLParams{Basis: badHammingBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.verifyMoufang2()
	if err == nil {
		t.Fatalf("Bad Hamming basis passed Moufang test.")
	}
}

// This is commented out because it is VERY SLOW

func TestGolayMoufang2(t *testing.T) {
	cl, err := NewCL(CLParams{Basis: GolayBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.verifyMoufang2()
	if err != nil {
		t.Fatalf("Golay basis failed verifyMoufang(): %s", err)
	}
}

func TestBadGolayNotMoufang(t *testing.T) {
	cl, err := NewCL(CLParams{Basis: badGolayBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.verifyMoufang()
	if err == nil {
		t.Fatalf("Bad Golay basis passed verifyMoufang(), expected it to fail.")
	}
}

func TestPrintHammingElems(t *testing.T) {
	cl, err := NewCL(CLParams{Basis: HammingBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	cl.PrintLoopElems()
}

func TestPrintHammingVectorSpace(t *testing.T) {
	cl, err := NewCL(CLParams{Basis: HammingBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	cl.PrintVectorSpace()
}

func TestHammingTheta(t *testing.T) {

	for i := 0; i < 1000; i++ {
		basis := HammingBasis
		cl, err := NewCL(CLParams{Basis: basis, Theta: RandomTheta})
		if err != nil {
			t.Fatalf("Failed to create CL: %s", err)
		}
		vs := cl.VectorSpace()

		// Verify that the cocycle is normalized, then check the axioms (S), (C),
		// (A) from [Gri86] p. 225

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
				t.Fatalf("Expected theta(v%d,v%d) to be %d, got %d", i-1, i-1, (BitWeight(vs[i])/4)%2, b)
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
}

// func TestGolayTheta(t *testing.T) {

// 	basis := GolayBasis
// 	cl, err := NewCL(CLParams{Basis: basis)
// 	if err != nil {
// 		t.Fatalf("Failed to create CL: %s", err)
// 	}
// 	vs := cl.VectorSpace()

// 	// Verify that the cocycle is normalized, then check the axioms (S), (C),
// 	// (A) from [Gri86] p. 225

// 	// theta(0,x) == theta(x,0) == 0 (normalized cocycle)
// 	for _, v := range vs {
// 		if b, _ := cl.ThetaByVec(v, 0); b != 0 {
// 			t.Fatalf("Theta not normalized at %x, 0", v)
// 		}
// 		if b, _ := cl.ThetaByVec(0, v); b != 0 {
// 			t.Fatalf("Theta not normalized at 0, %x", v)
// 		}
// 	}

// 	// (S) for all x, theta(x,x) === |x|/4 (all congruence is mod 2)
// 	for i := 0; i < len(vs); i++ {
// 		if b, _ := cl.ThetaByIdx(uint(i), uint(i)); b != (BitWeight(vs[i])/4)%2 {
// 			t.Fatalf("Expected theta(v%d,v%d) to be %d, got %d", i-1, i-1, (BitWeight(vs[i])/4)%2, b)
// 		}
// 	}

// 	// (C) for all x,y, theta(x,y) + theta(y,x) === |x&y|/2
// 	for i := 0; i < len(vs); i++ {
// 		for j := 0; j < len(vs); j++ {
// 			x, y := vs[i], vs[j]
// 			a, _ := cl.ThetaByIdx(uint(i), uint(j))
// 			b, _ := cl.ThetaByIdx(uint(j), uint(i))
// 			lhs := uint((a + b) % 2)
// 			rhs := (BitWeight(x&y) / 2) % 2
// 			if rhs != lhs {
// 				t.Fatalf("Expected theta(%x,%x) + theta(%x,%x) to be %d, got %d", x, y, y, x, rhs, lhs)
// 			}
// 		}
// 	}

// 	// (A) for all x,y,z theta(x,y) + theta(x^y,z) + theta(y,z) + theta(x,y^z) === |x&y&z|
// 	for i := 0; i < len(vs); i++ {
// 		for j := 0; j < len(vs); j++ {
// 			for k := 0; k < len(vs); k++ {
// 				x, y, z := vs[i], vs[j], vs[k]
// 				a, _ := cl.ThetaByVec(x, y)
// 				b, _ := cl.ThetaByVec(x^y, z)
// 				c, _ := cl.ThetaByVec(y, z)
// 				d, _ := cl.ThetaByVec(x, y^z)
// 				lhs := (a + b + c + d) % 2
// 				rhs := (BitWeight(x & y & z)) % 2
// 				if rhs != lhs {
// 					t.Errorf("(x,y - %x,%x): %d (x^y,z - %x,%x): %d (x,y^z - %x,%x): %d (y,z - %x,%x): %d\n",
// 						x, y, a,
// 						x^y, z, b,
// 						x, y^z, c,
// 						y, z, d)
// 					t.Fatalf("Error in triple identity for %x %x %x, Expected %d, got %d", x, y, z, rhs, lhs)
// 				}
// 			}
// 		}
// 	}
// }

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
	cl, err := NewCL(CLParams{Basis: HammingBasis})
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
	cl, err := NewCL(CLParams{Basis: GolayBasis})
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
	cl, err := NewCL(CLParams{Basis: HammingBasis})
	if err != nil {
		b.Fatalf("Failed to create CL: %s", err)
	}
	for n := 0; n < b.N; n++ {
		err = cl.verifyMoufang()
		if err != nil {
			b.Fatalf("%s", err)
		}
	}
}

func BenchmarkHammingMoufang2(b *testing.B) {
	cl, err := NewCL(CLParams{Basis: HammingBasis})
	if err != nil {
		b.Fatalf("Failed to create CL: %s", err)
	}
	for n := 0; n < b.N; n++ {
		err = cl.verifyMoufang2()
		if err != nil {
			b.Fatalf("%s", err)
		}
	}
}

func BenchmarkMul(b *testing.B) {
	cl, err := NewCL(CLParams{Basis: HammingBasis})
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
	cl, err := NewCL(CLParams{Basis: HammingBasis})
	if err != nil {
		b.Fatalf("Failed to create CL: %s", err)
	}
	for n := 0; n < b.N; n++ {
		cl.LoopElems()
	}
}

func BenchmarkLoopElemsGolay(b *testing.B) {
	cl, err := NewCL(CLParams{Basis: GolayBasis})
	if err != nil {
		b.Fatalf("Failed to create CL: %s", err)
	}
	for n := 0; n < b.N; n++ {
		cl.LoopElems()
	}
}
