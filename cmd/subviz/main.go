package main

import (
	"fmt"
	"github.com/bnagy/codeloops"
	"github.com/fogleman/gg"
	"log"
	"math"
)

const PIXEL = 16 // This should be 0 % 8

// don't bother drawing a label below a certain pixel size - the text just
// greeks
const MIN_LABEL_H = 10

// Likewise, we don't need a VAST label above the 4096x4096 cocycle
const MAX_LABEL_H = 50

// Do we want to label the elements in the multiplication table?
const CELL_LABELS = false

func main() {

	basis := codeloops.GolayBasis
	choose := 7

	subspaces := 0
	codeloops.SetCombinationsWithoutReplacement(uint(len(basis)), uint(choose), func(s []uint) {
		// such a lazy way to do this, but has the advantage that the number
		// here will definitely match the loop below
		subspaces++
	})
	xreps := int(math.Ceil(math.Sqrt(float64(subspaces))))
	yreps := xreps
	// but we only display the +/+ top left corner
	// since the others are copies (two sign flipped)
	imgInner := (1 << uint(choose)) * PIXEL // * 2 for full loop
	// How many actual loop elems are we drawing?
	pxSquare := (1 << uint(choose)) // * 2 for full loop

	// I am from the 'moar variables moar better' school.
	imgInnerf := float64(imgInner)
	borderW := PIXEL / 8
	if borderW == 0 {
		borderW = PIXEL
	}
	borderWf := float64(borderW)
	borderLenf := imgInnerf + 2*borderWf
	imgW := imgInner + 4*int(borderW) // black border, white border * 2 sides
	imgWf := float64(imgW)
	labelHf := math.Ceil(imgWf / 6) // arbitrarily pleasing ratio
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
	dc := gg.NewContext(panelW*xreps, panelH*yreps)
	dc.SetRGB(1, 1, 1)
	dc.Clear()
	err := dc.LoadFontFace("/Library/Fonts/Envy Code R.ttf", labelHf*3/8)
	if err != nil {
		log.Fatalf("Unable to load font Envy Code R: %s", err)
	}

	count := 0
	codeloops.SetCombinationsWithoutReplacement(uint(len(basis)), uint(choose), func(s []uint) {

		dc.SetRGB(0, 0, 0)
		tlXf := float64(count%xreps) * panelWf // top left X coordinate of the panel
		tlYf := float64(count/yreps) * panelHf
		count++

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
		thisBasis := []uint{}
		for _, idx := range s {
			thisBasis = append(thisBasis, codeloops.GolayBasis[idx])
		}
		cl, _ := codeloops.NewCL(codeloops.CLParams{Basis: thisBasis})

		if labelH > 0 {
			label := "Basis:"
			for _, idx := range s {
				label += fmt.Sprintf(" v%d", idx+1)
			}
			// Drawing two strings in a whitespace height of labelHf+borderWf,
			// so one is vertically centred at h/4, one at 3h/4
			dc.DrawStringAnchored(label, tlXf+borderWf, tlYf+(labelHf+borderWf)/4, 0, 0.5)
			dc.DrawStringAnchored(fmt.Sprintf("Seed: %s", cl.Seed), tlXf+borderWf, tlYf+(labelHf+borderWf)*3/4, 0, 0.5)
		}

		res := new(codeloops.CLElem)
		vm := cl.VectorIdxMap()
		loopElems := cl.LoopElems()
		negColor, posColor := "", ""
		if cl.IsAssoc() {
			negColor = "61a62f" // red and green
			posColor = "ae2013"
		} else if cl.IsMoufang() {
			negColor = "3b57e2" // blue and yellow
			posColor = "f8cc74"
		} else {
			// if this happens, it is either a bug or a mathematical
			// discovery.
			negColor = "6F824A" // puke green and rotten green
			posColor = "8dae4b"
		}
		for i := 0; i < int(pxSquare); i++ {
			for j := 0; j < int(pxSquare); j++ {
				_, e := cl.Mul(&loopElems[i], &loopElems[j], res)
				if e != nil {
					log.Fatal(e)
				}
				if res.Sign() > 0 {
					dc.SetHexColor(posColor)
				} else {
					dc.SetHexColor(negColor)
				}
				dc.DrawRectangle(
					tlXf+2*borderWf+float64(j*PIXEL),         // x
					tlYf+labelHf+2*borderWf+float64(i*PIXEL), // y
					float64(PIXEL),                           // w
					float64(PIXEL))                           // h
				dc.Fill()
				if CELL_LABELS {
					dc.SetRGB(0, 0, 0)
					dc.DrawStringAnchored(
						fmt.Sprintf("%x", vm[res.Vec()]),
						tlXf+2*borderWf+float64(j*PIXEL)+PIXEL/2,
						tlYf+labelHf+2*borderWf+float64(i*PIXEL)+PIXEL/2,
						0.5,
						0.5)
				}
			}
		}
	})

	dc.SavePNG("out.png")
}
