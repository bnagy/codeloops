package main

import (
	"github.com/bnagy/codeloops"
	"github.com/fogleman/gg"
	"log"
)

const PIXEL = 40 // This should be even!
// const alphaLen = 128

func main() {

	// Create the alpha space, and build a full theta from the corresponding
	// basis. We are building two spaces V and W such that any vector in the
	// Golay space (whose basis is 12 vectors) can be written as v + w for
	// some v in V, w in W. This particular basis also has the additional
	// property that the subspace V forms a group, so theta(v1,v2) is always
	// 0. See the partition code for details.

	// Which basis to use? Hopefully this is the only thing that needs
	// configuring, but the border might look ugly for very small bases.
	B := codeloops.GolayAwesumBasis

	basisLen := len(B)
	alphaLen := 1 << (uint(basisLen)/2 + 1)
	v := B[:basisLen/2]
	w := B[basisLen/2:]
	cl, _ := codeloops.NewCL(codeloops.CLParams{Basis: B})
	alpha := []uint{}
	vsV := codeloops.VectorSpace(v)
	vsW := codeloops.VectorSpace(w) // [1:] // we don't want the zero vector at the start
	alpha = append(vsV, vsW...)

	//alphaLen = alphaLen / 2
	// Set up image size
	imgInner := alphaLen * PIXEL
	imgInnerf := float64(imgInner)
	borderW := imgInner / 100
	if borderW == 0 {
		borderW = PIXEL
	}
	borderWf := float64(borderW)
	borderLenf := imgInnerf + 2*borderWf
	imgW := imgInner + 4*int(borderW) // black border, white border * 2 sides
	imgH := imgInner + 4*int(borderW) // + int(labelHf)

	// Set up drawing context
	dc := gg.NewContext(imgW, imgH)
	dc.SetRGB(1, 1, 1)
	dc.Clear()
	dc.SetRGB(0, 0, 0)
	tlXf := float64(0)
	tlYf := float64(0)

	// draw the top border
	dc.DrawRectangle(tlXf+borderWf, tlYf+borderWf, borderLenf, borderWf)
	dc.Fill()
	// left...
	dc.DrawRectangle(tlXf+borderWf, tlYf+borderWf, borderWf, borderLenf)
	dc.Fill()
	// bottom...
	dc.DrawRectangle(tlXf+borderWf, tlYf+borderWf*2+imgInnerf, borderLenf, borderWf)
	dc.Fill()
	// right.
	dc.DrawRectangle(tlXf+borderWf*2+imgInnerf, tlYf+borderWf, borderWf, borderLenf)
	dc.Fill()

	// draw a black rectangle for all theta(x,y) == 1
	for i := 0; i < alphaLen; i++ {
		for j := 0; j < alphaLen; j++ {
			res, e := cl.ThetaByVec(alpha[i], alpha[j])
			if e != nil {
				log.Fatal(e)
			}
			if res > 0 {
				dc.DrawRectangle(
					tlXf+2*borderWf+float64(j*PIXEL), // x
					tlYf+2*borderWf+float64(i*PIXEL), // y
					float64(PIXEL),                   // w
					float64(PIXEL))                   // h
				dc.Fill()
			}
		}
	}

	dc.SavePNG("alpha.png")
}
