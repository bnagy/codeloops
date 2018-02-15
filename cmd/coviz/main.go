package main

import (
	"fmt"
	"github.com/bnagy/codeloops"
	"github.com/fogleman/gg"
	"log"
	"math"
)

const PIXEL = 8 // This should be even!
const XREPS = 40
const YREPS = 20

// don't bother drawing a label below a certain pixel size - the text just
// greeks
const MIN_LABEL_H = 10

// Likewise, we don't need a VAST label above the 4096x4096 cocycle
const MAX_LABEL_H = 50

func main() {

	basis := codeloops.HammingBasis
	cl, err := codeloops.NewCL(codeloops.CLParams{Basis: basis})
	if err != nil {
		log.Fatalf("Failed to create CL: %s", err)
	}
	vs := cl.VectorSpace()

	imgInner := (1 << uint(len(basis)) * PIXEL)
	imgInnerf := float64(imgInner)
	borderW := PIXEL / 2
	if borderW == 0 {
		borderW = PIXEL
	}
	borderWf := float64(borderW)
	borderLenf := imgInnerf + 2*borderWf
	imgW := imgInner + 4*int(borderW) // black border, white border * 2 sides
	imgWf := float64(imgW)
	labelHf := math.Ceil(imgWf / 12) // arbitrarily pleasing ratio
	if float64(imgW/12) > float64(MAX_LABEL_H) {
		labelHf = float64(MAX_LABEL_H)
	}
	if float64(imgW/12) < float64(MIN_LABEL_H) {
		labelHf = 0
	}
	labelH := int(labelHf) // at small PIXELS this will be 0
	imgH := imgInner + 4*int(borderW) + int(labelHf)
	// complete WxH of one panel (image + label + borders + gutter)
	panelH := imgH
	panelW := imgW
	panelHf, panelWf := float64(panelH), float64(panelW)

	dc := gg.NewContext(panelW*XREPS, panelH*YREPS)
	dc.SetRGB(1, 1, 1)
	dc.Clear()
	err = dc.LoadFontFace("/Library/Fonts/Envy Code R.ttf", labelHf)
	if err != nil {
		log.Fatalf("Unable to load font Envy Code R: %s", err)
	}

	i := uint(0x2500)
	for yrep := float64(0); yrep < YREPS; yrep++ {

		for xrep := float64(0); xrep < XREPS; xrep++ {

			cl, _ := codeloops.NewCL(codeloops.CLParams{Basis: basis, Theta: i})
			i += 8 // adding 1 is a nop because the first path choice is always 0
			dc.SetRGB(0, 0, 0)
			tlXf := xrep * panelWf // top left point of the panel
			tlYf := yrep * panelHf

			// top border
			dc.DrawRectangle(tlXf+borderWf, tlYf+labelHf+borderWf, borderLenf, borderWf)
			dc.Fill()
			// left...
			dc.DrawRectangle(tlXf+borderWf, tlYf+labelHf+borderWf, borderWf, borderLenf)
			dc.Fill()
			// bottom...
			dc.DrawRectangle(tlXf+borderWf, tlYf+labelHf+borderWf*2+imgInnerf, borderLenf, borderWf)
			dc.Fill()
			// right.
			dc.DrawRectangle(tlXf+borderWf*2+imgInnerf, tlYf+labelHf+borderWf, borderWf, borderLenf)
			dc.Fill()

			// draw the label
			if labelH > 0 {
				label := fmt.Sprintf("Theta: %s", cl.ThetaPath)
				// Center horizontally in the space made of white border + label gutter
				dc.DrawStringAnchored(label, tlXf+borderWf, tlYf+labelHf/2+borderWf/2, 0, 0.5)
			}
			// draw a black rectangle for all theta(x,y) == 1
			for i, _ := range vs {
				for j, _ := range vs {
					res, e := cl.ThetaByIdx(uint(i), uint(j))
					if e != nil {
						log.Fatal(e)
					}
					if res > 0 {
						dc.DrawRectangle(
							tlXf+2*borderWf+float64(j*PIXEL),         // x
							tlYf+labelHf+2*borderWf+float64(i*PIXEL), // y
							float64(PIXEL),                           // w
							float64(PIXEL))                           // h
						dc.Fill()
					}
				}
			}
		}
	}

	dc.SavePNG("out.png")
}
