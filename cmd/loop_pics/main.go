package main

import (
	"github.com/bnagy/codeloops"
	"github.com/fogleman/gg"
	"log"
)

const PIXEL = 2    // This should be even!
const TRUNCATE = 0 // at this many pixels

func SpecialVectorSpace(b6, b5, b1 []uint) []uint {
	// We assume that b6 and b5 generate subspaces with only 0 in common, and
	// that b1 is not in vs6 or vs5.

	final := []uint{}
	vsV := codeloops.VectorSpace(b6)
	bW := append(b5, b1...)
	vsW := codeloops.VectorSpace(bW)
	final = append(final, vsV...)

	for _, v := range vsV {
		for _, w := range vsW[1:] {
			final = append(final, v^w)
		}
	}

	if len(final) != 1<<12 {
		log.Fatalf("Wrong number of vectors constructing special vector space ordering: %d", len(final))
	}
	return final
}

func main() {

	B := codeloops.GolayAwesumBasis

	var pxSquare uint
	var imgInner int
	if TRUNCATE > 0 {
		imgInner = (TRUNCATE * PIXEL)
		pxSquare = TRUNCATE
	} else {
		imgInner = (1 << uint(len(B)) * PIXEL)
		pxSquare = 1 << uint(len(B))
	}

	// I am from the 'moar variables moar better' school.
	imgInnerf := float64(imgInner)
	borderW := PIXEL / 2
	if borderW == 0 {
		borderW = PIXEL
	}
	borderWf := float64(borderW)
	borderLenf := imgInnerf + 2*borderWf
	imgW := imgInner + 4*int(borderW) // black border, white border * 2 sides
	imgH := imgInner + 4*int(borderW)

	dc := gg.NewContext(imgW, imgH)
	dc.SetRGB(1, 1, 1)
	dc.Clear()

	cl, _ := codeloops.NewCL(codeloops.CLParams{Basis: B})
	// This will only work for bases that have the special partitioning
	// property (see cmd/partition). Otherwise, use the line below.
	vs := SpecialVectorSpace(B[:6], B[6:11], B[11:12])
	// vs := codeloops.VectorSpace(B)

	dc.SetRGB(0, 0, 0)

	// top border
	dc.DrawRectangle(borderWf, borderWf, borderLenf, borderWf)
	dc.Fill()
	// left...
	dc.DrawRectangle(borderWf, borderWf, borderWf, borderLenf)
	dc.Fill()
	// bottom...
	dc.DrawRectangle(borderWf, borderWf*2+imgInnerf, borderLenf, borderWf)
	dc.Fill()
	// right.
	dc.DrawRectangle(borderWf*2+imgInnerf, borderWf, borderWf, borderLenf)
	dc.Fill()

	// draw a black rectangle for all theta(x,y) == 1
	for i := 0; i < int(pxSquare); i++ {
		for j := 0; j < int(pxSquare); j++ {
			res, e := cl.ThetaByVec(vs[i], vs[j])
			if e != nil {
				log.Fatal(e)
			}
			if res > 0 {
				dc.DrawRectangle(
					2*borderWf+float64(j*PIXEL), // x
					2*borderWf+float64(i*PIXEL), // y
					float64(PIXEL),              // w
					float64(PIXEL))              // h
				dc.Fill()
			}
		}
	}

	dc.SavePNG("out.png")
}
