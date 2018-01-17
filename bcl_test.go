package codeloops

import (
	"runtime"
	"strconv"
	"strings"
	"sync"
	"testing"
)

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
	// populate an array with all elements of the extension
	elems := cl.LoopElems()
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
			t.Fatalf("Failed Moufang Identity for %s %s %s", x.String(), y.String(), z.String())
		}

	}

}

func BenchmarkHammingMoufang(b *testing.B) {
	for n := 0; n < b.N; n++ {
		cl, err := NewCL(HammingBasis)
		if err != nil {
			b.Fatalf("Failed to create CL: %s", err)
		}
		// populate an array with all elements of the extension
		elems := cl.LoopElems()
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

			// LHS
			z.Mul2(y, zy)
			x.Mul2(zy, x_zy)
			z.Mul2(x_zy, z__x_zy)
			//RHS
			z.Mul2(x, zx)
			zx.Mul2(z, zx_z)
			zx_z.Mul2(y, zx_z__y)

			if (z__x_zy.vec<<1)&z__x_zy.sgn != (zx_z__y.vec<<1)&zx_z__y.sgn {
				b.Fatalf("Failed Moufang Identity for %s %s %s", x.String(), y.String(), z.String())
			}
		}

	}
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
	zy, x_zy, z__x_zy := new(CLElem), new(CLElem), new(CLElem)
	zx, zx_z, zx_z__y := new(CLElem), new(CLElem), new(CLElem)

	// Iterate! Treat a bit string as 3 concatenated element indices, and then
	// use those to select the appropriate real element out of the work.Elems
	// slice.

	// Create masks to pull the element indicies from the counter
	elemSz := int(work.BaseCL.basisLen + 1)
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

		// verify Moufang identity:
		// z(x(zy)) = ((zx)z)y
		// These identifiers are not good Go style, but they capture the
		// associativity which is more important right now.
		// LHS:
		z.Mul2(y, zy)
		x.Mul2(zy, x_zy)
		z.Mul2(x_zy, z__x_zy)
		// RHS:
		z.Mul2(x, zx)
		zx.Mul2(z, zx_z)
		zx_z.Mul2(y, zx_z__y)
		// combine vec and sgn so we can compare both
		if (z__x_zy.vec<<1)&z__x_zy.sgn != (zx_z__y.vec<<1)&zx_z__y.sgn {
			work.Failures <- []*CLElem{x, y, z}
		}
	}
}

func TestGolayMoufangParallel(t *testing.T) {

	cl, err := NewCL(GolayBasis)
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}

	// populate an array with all elements of the extension
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
		go moufangParallelWorker(wu)
	}

	// Wait for the work to finish, or bail immediately if one of the workers
	// reports a failed triple.
loop:
	for {
		select {
		case fail := <-failures:
			t.Fatalf("Failed Moufang Identity for %s %s %s", fail[0].String(), fail[1].String(), fail[2].String())
		case <-done:
			// will fire when the done chan is closed
			break loop
		}
	}

}

func BenchmarkHammingMoufangParallel(b *testing.B) {

	for n := 0; n < b.N; n++ {
		cl, err := NewCL(HammingBasis)
		if err != nil {
			b.Fatalf("Failed to create CL: %s", err)
		}

		// populate an array with all elements of the extension
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
			go moufangParallelWorker(wu)
		}

		// Wait for the work to finish, or bail immediately if one of the workers
		// reports a failed triple.
	loop:
		for {
			select {
			case fail := <-failures:
				b.Fatalf("Failed Moufang Identity for %s %s %s", fail[0].String(), fail[1].String(), fail[2].String())
			case <-done:
				// will fire when the done chan is closed
				break loop
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
