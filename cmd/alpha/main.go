package main

import (
	"github.com/bnagy/codeloops"
	"github.com/fogleman/gg"
	"log"
)

const PIXEL = 4 // This should be even!
const ALPHA = 127

func main() {

	// Create the alpha space, and build a full theta from the corresponding
	// basis. We are building two spaces V and W such that any vector in the
	// Golay space (whose basis is 12 vectors) can be written as v + w for
	// some v in V, w in W. This particular basis also has the additional
	// property that the subspace V forms a group, so theta(v1,v2) is always
	// 0. See the partition code for details.
	v := codeloops.GolaySplitBasis[:6]
	w := codeloops.GolaySplitBasis[6:]
	cl, _ := codeloops.NewCL(codeloops.CLParams{Basis: codeloops.GolaySplitBasis})
	alpha := []uint{}
	vsV := codeloops.VectorSpace(v)
	vsW := codeloops.VectorSpace(w)[1:] // we don't want the zero vector at the start
	alpha = append(vsV, vsW...)
	if len(alpha) != ALPHA {
		log.Fatalf("Bad length for alpha space %d, expected %d", len(alpha), ALPHA)
	}

	// Set up image size
	imgInner := ALPHA * PIXEL
	imgInnerf := float64(imgInner)
	borderW := PIXEL / 2
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
	for i := 0; i < ALPHA; i++ {
		for j := 0; j < ALPHA; j++ {
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
