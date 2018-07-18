package codeloops

import (
	//"fmt"
	"testing"
)

func TestAlphaIdxMatchGolay(t *testing.T) {
	cl, err := NewCL(CLParams{Basis: GolaySplitBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	for i := uint(0); i < cl.alphaSz; i++ {
		for j := uint(0); j < cl.alphaSz; j++ {
			alphaRes, err := cl.AlphaByIdx(i, j)

			if err != nil {
				t.Fatalf("Failed AlphaByIdx: %s", err)
			}
			thetaRes, err := cl.ThetaByVec(cl.vsAlpha[i], cl.vsAlpha[j])
			if err != nil {
				t.Fatalf("Failed ThetaByVec: %s", err)
			}
			if alphaRes != thetaRes {
				t.Fatalf("Mismatch at (%d, %d): Alpha: %d, Theta: %d", i, j, alphaRes, thetaRes)
			}
		}
	}
}

func TestAlphaVecMatchGolay(t *testing.T) {
	cl, err := NewCL(CLParams{Basis: GolaySplitBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	for _, v := range cl.vsAlpha {
		for _, w := range cl.vsAlpha {
			alphaRes, err := cl.AlphaByVec(v, w)

			if err != nil {
				t.Fatalf("Failed AlphaByVec: %s", err)
			}
			thetaRes, err := cl.ThetaByVec(v, w)
			if err != nil {
				t.Fatalf("Failed ThetaByVec: %s", err)
			}
			if alphaRes != thetaRes {
				t.Fatalf("Mismatch with 0x%x, 0x%x: Alpha: %d, Theta: %d", v, w, alphaRes, thetaRes)
			}
		}
	}
}

func TestAlphaIdxMatchHamming(t *testing.T) {
	cl, err := NewCL(CLParams{Basis: HammingBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	for i := uint(0); i < cl.alphaSz; i++ {
		for j := uint(0); j < cl.alphaSz; j++ {
			alphaRes, err := cl.AlphaByIdx(i, j)

			if err != nil {
				t.Fatalf("Failed AlphaByIdx: %s", err)
			}
			thetaRes, err := cl.ThetaByVec(cl.vsAlpha[i], cl.vsAlpha[j])
			if err != nil {
				t.Fatalf("Failed ThetaByVec: %s", err)
			}
			if alphaRes != thetaRes {
				t.Fatalf("Mismatch at (%d, %d): Alpha: %d, Theta: %d", i, j, alphaRes, thetaRes)
			}
		}
	}
}

func TestAlphaVecMatchHamming(t *testing.T) {
	cl, err := NewCL(CLParams{Basis: HammingBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	for _, v := range cl.vsAlpha {
		for _, w := range cl.vsAlpha {
			alphaRes, err := cl.AlphaByVec(v, w)

			if err != nil {
				t.Fatalf("Failed AlphaByVec: %s", err)
			}
			thetaRes, err := cl.ThetaByVec(v, w)
			if err != nil {
				t.Fatalf("Failed ThetaByVec: %s", err)
			}
			if alphaRes != thetaRes {
				t.Fatalf("Mismatch with 0x%x, 0x%x: Alpha: %d, Theta: %d", v, w, alphaRes, thetaRes)
			}
		}
	}
}

func TestAlphaHamming(t *testing.T) {
	cl, err := NewCL(CLParams{Basis: HammingBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	for _, v := range cl.vs {
		for _, w := range cl.vs {
			alphaRes, err := cl.ThetaAlphaByVec(v, w)
			if err != nil {
				t.Fatalf("Failed in ThetaAlpha: %s", err)
			}
			thetaRes, err := cl.ThetaByVec(v, w)
			if err != nil {
				t.Fatalf("Failed in Theta: %s", err)
			}
			if alphaRes != thetaRes {
				t.Fatalf("Mismatch with 0x%x, 0x%x: Alpha: %d, Theta: %d", v, w, alphaRes, thetaRes)
			}
		}
	}
}

func TestAlphaGolay(t *testing.T) {
	cl, err := NewCL(CLParams{Basis: GolaySplitBasis})
	if err != nil {
		t.Fatalf("Failed to create CL: %s", err)
	}
	for _, v := range cl.vs {
		for _, w := range cl.vs {
			alphaRes, err := cl.ThetaAlphaByVec(v, w)
			if err != nil {
				t.Fatalf("Failed in ThetaAlpha: %s", err)
			}
			thetaRes, err := cl.ThetaByVec(v, w)
			if err != nil {
				t.Fatalf("Failed in Theta: %s", err)
			}
			if alphaRes != thetaRes {
				t.Fatalf("Mismatch with 0x%x, 0x%x: Alpha: %d, Theta: %d", v, w, alphaRes, thetaRes)
			}
		}
	}
}

func BenchmarkAlphaGolay(b *testing.B) {
	cl, err := NewCL(CLParams{Basis: GolayBasis})
	if err != nil {
		b.Fatalf("Failed to create CL: %s", err)
	}
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		v := cl.vs[b.N%len(cl.vs)]
		w := cl.vs[len(cl.vs)-b.N%len(cl.vs)-1]
		_, err = cl.ThetaAlphaByVec(v, w)
		if err != nil {
			b.Fatalf("Failed in ThetaAlpha: %s", err)
		}

	}
}

func BenchmarkAlphaGolayFast(b *testing.B) {
	cl, err := NewCL(CLParams{Basis: GolayBasis})
	if err != nil {
		b.Fatalf("Failed to create CL: %s", err)
	}
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		v := cl.vs[b.N%len(cl.vs)]
		w := cl.vs[len(cl.vs)-b.N%len(cl.vs)-1]
		cl.thetaAlphaByVecFast(v, w)
	}
}

func BenchmarkThetaGolay(b *testing.B) {
	cl, err := NewCL(CLParams{Basis: GolayBasis})
	if err != nil {
		b.Fatalf("Failed to create CL: %s", err)
	}
	b.ResetTimer()
	for n := 0; n < b.N; n++ {
		v := cl.vs[b.N%len(cl.vs)]
		w := cl.vs[len(cl.vs)-b.N%len(cl.vs)-1]
		_, err = cl.ThetaByVec(v, w)
		if err != nil {
			b.Fatalf("Failed in Theta: %s", err)
		}
	}
}
