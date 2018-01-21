package codeloops

import (
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
		if bitWeight(x) != expect[idx] {
			t.Errorf("TestBitWeight: Weight of 0x%x should be %d, got %d", x, expect[idx], bitWeight(x))
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

func TestVerifyGolayMoufangGood(t *testing.T) {
	cl, err := NewCL(GolayBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	err = cl.VerifyMoufang()
	if err != nil {
		t.Fatalf("Golay basis failed VerifyMoufang(): %s", err)
	}
}

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

func BenchmarkBitWeight(b *testing.B) {
	// dummy var stops the compiler from eliminating dead code in the loop
	dummy := uint(0)
	for i := 0; i < b.N; i++ {
		dummy = bitWeight(0xdeadbeef)
	}
}

// This benches much slower, despite being implemented as a single asm
// instruction? Maybe the trampoline process costs a nanosecond or so...
func BenchmarkPopCntq(b *testing.B) {
	dummy := uint(0)
	for i := 0; i < b.N; i++ {
		dummy = popCntq(0xdeadbeef)
	}
}

func BenchmarkSigma(b *testing.B) {
	dummy := uint(0)
	for i := 0; i < b.N; i++ {
		dummy = sigma(0xdeadbeef, 0xfacef00d)
	}
}

func BenchmarkVerifyBasisHamming(b *testing.B) {
	for n := 0; n < b.N; n++ {
		cl, err := NewCL(HammingBasis)
		if err != nil {
			b.Fatalf("Failed to create CL: %s", err)
		}
		err = cl.VerifyBasis()
		if err != nil {
			b.Fatalf("Hamming basis failed VerifyBasis(), expected it to pass.")
		}
	}
}

func BenchmarkVerifyBasisGolay(b *testing.B) {
	for n := 0; n < b.N; n++ {
		cl, err := NewCL(GolayBasis)
		if err != nil {
			b.Fatalf("Failed to create CL: %s", err)
		}
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
