package codeloops

import (
	"fmt"
	"strconv"
	"strings"
	"testing"
)

func nextProduct(a []CLElem, k int) func() []CLElem {
	// ref: https://stackoverflow.com/questions/23412146/create-cartesian-product-in-go
	p := make([]CLElem, k)
	x := make([]int, len(p))
	return func() []CLElem {
		p = p[:len(x)]
		for i, xi := range x {
			p[i] = a[xi]
		}
		for i := len(x) - 1; i >= 0; i-- {
			x[i]++
			if x[i] < len(a) {
				break
			}
			x[i] = 0
			if i <= 0 {
				x = x[0:0]
				break
			}
		}
		return p
	}
}

func createIterator(cl *CL, k int) (np func() []CLElem, e error) {

	// populate an array with all elements of the extension
	elems := []CLElem{}
	for i := uint(0); i < cl.size; i++ {
		ipos, err := cl.NewCLElem(i, Pos)
		if err != nil {
			e = fmt.Errorf("Error in NewCLElem(%d): %s", ipos, err)
			return
		}
		elems = append(elems, *ipos)
		ineg, err := cl.NewCLElem(i, Neg)
		if err != nil {
			e = fmt.Errorf("Error in NewCLElem(%d): %s", ineg, err)
			return
		}
		elems = append(elems, *ineg)
	}

	// return the n-tuple iterator function
	np = nextProduct(elems, k)
	return
}

func TestBitWeight(t *testing.T) {
	tests := []uint{0x0, 0x1, 0x87, 0xff, 0xffffffff, 0xdeadbeef, 0xffffffffdeadbeef}
	expect := []int{0, 1, 4, 8, 32, 24, 56}
	for idx, x := range tests {
		if bitWeight(x) != expect[idx] {
			t.Errorf("TestBitWeight: Weight of 0x%x should be %d, got %d", x, expect[idx], bitWeight(x))
		}
	}
}

func TestParity(t *testing.T) {
	tests := []uint{0x0, 0x1, 0x87, 0xff, 0xffffffff, 0xdeadbeef, 0xefffffffdeadbeef}
	expect := []int{0, 1, 0, 0, 0, 0, 1}
	for idx, x := range tests {
		if parity(x) != expect[idx] {
			t.Errorf("TestParity: Parity of 0x%x should be %d, got %d", x, expect[idx], bitWeight(x))
		}
	}
}

func TestHammingMoufang(t *testing.T) {
	cl, err := NewCL(HammingBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}

	// Lazily deal with all triples from the extension set
	iter, err := createIterator(cl, 3)
	if err != nil {
		t.Fatalf("Error in TestHammingMoufang creating iterator: %s", err)
	}
	errors := 0
	for {
		product := iter()
		if len(product) == 0 {
			break
		}
		// verify Moufang identity:
		// z(x(zy)) = ((zx)z)y

		x, y, z := &product[0], &product[1], &product[2]

		// These identifiers are not good Go style, but they capture the associativity
		// which is more important right now...
		// LHS 1
		zy := z.Mul(y)
		x_zy := x.Mul(zy)
		z__x_zy := z.Mul(x_zy)
		//RHS 1
		zx := z.Mul(x)
		zx_z := zx.Mul(z)
		zx_z__y := zx_z.Mul(y)
		if z__x_zy.String() != zx_z__y.String() {
			t.Fatalf("Failed Moufang Identity for %s %s %s", x.String(), y.String(), z.String())
			errors++
		}

	}

}

func BenchmarkHammingMoufang(b *testing.B) {
	for n := 0; n < b.N; n++ {
		cl, err := NewCL(HammingBasis)
		if err != nil {
			b.Fatalf("Failed to create CL: %s", err)
		}

		// Lazily deal with all triples from the extension set
		iter, err := createIterator(cl, 3)
		if err != nil {
			b.Fatalf("Error in TestHammingMoufang creating iterator: %s", err)
		}
		errors := 0
		for {
			product := iter()
			if len(product) == 0 {
				break
			}
			// verify Moufang identity:
			// z(x(zy)) = ((zx)z)y

			x, y, z := &product[0], &product[1], &product[2]

			// These identifiers are not good Go style, but they capture the associativity
			// which is more important right now...
			// LHS 1
			zy := z.Mul(y)
			x_zy := x.Mul(zy)
			z__x_zy := z.Mul(x_zy)
			//RHS 1
			zx := z.Mul(x)
			zx_z := zx.Mul(z)
			zx_z__y := zx_z.Mul(y)
			if z__x_zy.String() != zx_z__y.String() {
				b.Fatalf("Failed Moufang Identity for %s %s %s", x.String(), y.String(), z.String())
				errors++
			}

		}
	}

}

func TestHammingMoufang2(t *testing.T) {
	cl, err := NewCL(HammingBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}

	// Lazily deal with all triples from the extension set
	iter, err := createIterator(cl, 3)
	if err != nil {
		t.Fatalf("Error in TestHammingMoufang creating iterator: %s", err)
	}
	errors := 0
	zy, x_zy, z__x_zy := new(CLElem), new(CLElem), new(CLElem)
	zx, zx_z, zx_z__y := new(CLElem), new(CLElem), new(CLElem)
	for {
		product := iter()
		if len(product) == 0 {
			break
		}
		// verify Moufang identity:
		// z(x(zy)) = ((zx)z)y

		x, y, z := &product[0], &product[1], &product[2]

		// These identifiers are not good Go style, but they capture the associativity
		// which is more important right now...
		// LHS 1
		z.Mul2(y, zy)
		x.Mul2(zy, x_zy)
		z.Mul2(x_zy, z__x_zy)
		//RHS 1
		z.Mul2(x, zx)
		zx.Mul2(z, zx_z)
		zx_z.Mul2(y, zx_z__y)
		if (z__x_zy.vec<<1)&z__x_zy.sgn != (zx_z__y.vec<<1)&zx_z__y.sgn {
			t.Fatalf("Failed Moufang Identity for %s %s %s", x.String(), y.String(), z.String())
			errors++
		}

	}
}

func BenchmarkHammingMoufang2(b *testing.B) {
	for n := 0; n < b.N; n++ {
		cl, err := NewCL(HammingBasis)
		if err != nil {
			b.Fatalf("Failed to create CL: %s", err)
		}

		// Lazily deal with all triples from the extension set
		iter, err := createIterator(cl, 3)
		if err != nil {
			b.Fatalf("Error in TestHammingMoufang creating iterator: %s", err)
		}
		x, y, z := new(CLElem), new(CLElem), new(CLElem)
		zy, x_zy, z__x_zy := new(CLElem), new(CLElem), new(CLElem)
		zx, zx_z, zx_z__y := new(CLElem), new(CLElem), new(CLElem)
		product := []CLElem{}
		for {
			product = iter()
			if len(product) == 0 {
				break
			}
			// verify Moufang identity:
			// z(x(zy)) = ((zx)z)y

			x, y, z = &product[0], &product[1], &product[2]

			// These identifiers are not good Go style, but they capture the associativity
			// which is more important right now...
			// LHS 1
			z.Mul2(y, zy)
			x.Mul2(zy, x_zy)
			z.Mul2(x_zy, z__x_zy)
			//RHS 1
			z.Mul2(x, zx)
			zx.Mul2(z, zx_z)
			zx_z.Mul2(y, zx_z__y)
			if (z__x_zy.vec<<1)&z__x_zy.sgn != (zx_z__y.vec<<1)&zx_z__y.sgn {
				b.Fatalf("Failed Moufang Identity for %s %s %s", x.String(), y.String(), z.String())
			}

		}
	}

}

func BenchmarkHammingMoufang3(b *testing.B) {
	for n := 0; n < b.N; n++ {
		cl, err := NewCL(HammingBasis)
		if err != nil {
			b.Fatalf("Failed to create CL: %s", err)
		}
		// populate an array with all elements of the extension
		elems := []CLElem{}
		for i := uint(0); i < cl.size; i++ {
			ipos, err := cl.NewCLElem(i, Pos)
			if err != nil {
				b.Fatalf("Error in NewCLElem(%d): %s", ipos, err)
			}
			elems = append(elems, *ipos)
			ineg, err := cl.NewCLElem(i, Neg)
			if err != nil {
				b.Fatalf("Error in NewCLElem(%d): %s", ineg, err)
			}
			elems = append(elems, *ineg)
		}
		// Preallocate the memory for the multiplication elems
		x, y, z := new(CLElem), new(CLElem), new(CLElem)
		zy, x_zy, z__x_zy := new(CLElem), new(CLElem), new(CLElem)
		zx, zx_z, zx_z__y := new(CLElem), new(CLElem), new(CLElem)

		// Iterate! Treat a bit string as 3 concatenated elems

		// Create masks to pull the elems from the counter
		elemSz := int(cl.basisLen + 1)
		zeroChunk := strings.Repeat("0", elemSz)
		oneChunk := strings.Repeat("1", elemSz)
		yMaskStr := zeroChunk + oneChunk + zeroChunk
		zMaskStr := zeroChunk + zeroChunk + oneChunk
		yMask, _ := strconv.ParseUint(yMaskStr, 2, 0)
		zMask, _ := strconv.ParseUint(zMaskStr, 2, 0)

		for combinedElems := uint(0); combinedElems < 1<<uint(elemSz*3); combinedElems++ {

			// verify Moufang identity:
			// z(x(zy)) = ((zx)z)y
			x = &elems[combinedElems>>10]              // top 5 bits as an index
			y = &elems[(combinedElems&uint(yMask))>>5] // middle 5 bits ...
			z = &elems[combinedElems&uint(zMask)]      // etc

			// These identifiers are not good Go style, but they capture the associativity
			// which is more important right now...
			// LHS 1
			z.Mul2(y, zy)
			x.Mul2(zy, x_zy)
			z.Mul2(x_zy, z__x_zy)
			//RHS 1
			z.Mul2(x, zx)
			zx.Mul2(z, zx_z)
			zx_z.Mul2(y, zx_z__y)
			if (z__x_zy.vec<<1)&z__x_zy.sgn != (zx_z__y.vec<<1)&zx_z__y.sgn {
				b.Fatalf("Failed Moufang Identity for %s %s %s", x.String(), y.String(), z.String())
			}

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
	for n := 0; n < b.N; n++ {
		x.Mul(y)
	}

}

func BenchmarkMul2(b *testing.B) {
	cl, err := NewCL(HammingBasis)
	if err != nil {
		b.Fatalf("Failed to create CL: %s", err)
	}
	x, _ := cl.NewCLElem(5, Pos)
	y, _ := cl.NewCLElem(7, Neg)
	z, _ := cl.NewCLElem(0, Pos)

	for n := 0; n < b.N; n++ {
		x.Mul2(y, z)
	}

}

func TestGolayMoufang(t *testing.T) {
	cl, err := NewCL(GolayBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	// populate an array with all elements of the extension
	elems := []CLElem{}
	for i := uint(0); i < cl.size; i++ {
		ipos, err := cl.NewCLElem(i, Pos)
		if err != nil {
			t.Fatalf("Error in NewCLElem(%d): %s", ipos, err)
		}
		elems = append(elems, *ipos)
		ineg, err := cl.NewCLElem(i, Neg)
		if err != nil {
			t.Fatalf("Error in NewCLElem(%d): %s", ineg, err)
		}
		elems = append(elems, *ineg)
	}
	// Preallocate the memory for the multiplication elems
	x, y, z := new(CLElem), new(CLElem), new(CLElem)
	zy, x_zy, z__x_zy := new(CLElem), new(CLElem), new(CLElem)
	zx, zx_z, zx_z__y := new(CLElem), new(CLElem), new(CLElem)

	// Iterate! Treat a bit string as 3 concatenated elems

	// Create masks to pull the elems from the counter
	elemSz := uint(cl.basisLen + 1)
	zeroChunk := strings.Repeat("0", int(elemSz))
	oneChunk := strings.Repeat("1", int(elemSz))
	yMaskStr := zeroChunk + oneChunk + zeroChunk
	zMaskStr := zeroChunk + zeroChunk + oneChunk
	yMask, _ := strconv.ParseUint(yMaskStr, 2, 0)
	zMask, _ := strconv.ParseUint(zMaskStr, 2, 0)

	for combinedElems := uint(0); combinedElems < 1<<uint(elemSz*3); combinedElems++ {

		// verify Moufang identity:
		// z(x(zy)) = ((zx)z)y

		//fmt.Printf("elems: %d\n", len(elems))
		x = &elems[combinedElems>>(elemSz*2)]           // top elemSz bits as an index
		y = &elems[(combinedElems&uint(yMask))>>elemSz] // middle elemSz bits ...
		z = &elems[combinedElems&uint(zMask)]           // etc

		// These identifiers are not good Go style, but they capture the associativity
		// which is more important right now...
		// LHS 1
		z.Mul2(y, zy)
		x.Mul2(zy, x_zy)
		z.Mul2(x_zy, z__x_zy)
		//RHS 1
		z.Mul2(x, zx)
		zx.Mul2(z, zx_z)
		zx_z.Mul2(y, zx_z__y)
		if (z__x_zy.vec<<1)&z__x_zy.sgn != (zx_z__y.vec<<1)&zx_z__y.sgn {
			t.Fatalf("Failed Moufang Identity for %s %s %s", x.String(), y.String(), z.String())
		}

	}

}
