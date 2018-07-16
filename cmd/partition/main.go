package main

import (
	"fmt"
	"github.com/bnagy/codeloops"
	"log"
)

func AddElems(to []int, from []uint, n int, verify func([]int) bool) (bool, []int) {
	for i := to[len(to)-1] + 1; i < len(from); i++ {
		candidate := append(to, i)
		if !verify(candidate) {
			// doesn't meet the criteria
			continue
		}

		// OK, we added this candidate and the set meets the criteria.
		if n == 1 {
			// this was the last one. Stop recursion.
			return true, candidate
		}
		// Still need n-1 more.
		if ok, set := AddElems(candidate, from, n-1, verify); ok {
			return true, set
		}

		// the candidate met the criteria, but we couldn't add enough extras.
		// Keep trying...
	}
	// did the whole loop and didn't find any candidates
	return false, nil
}

func main() {

	fullBasis := codeloops.GolayBasis
	fullVectorSpace := codeloops.VectorSpace(fullBasis)
	b6, b5 := []uint{}, []uint{}
	vs6, vs5 := []uint{}, []uint{}
	for i := 0; i < len(fullVectorSpace); i++ {
		s := []int{i}
		ok, set := AddElems(s, fullVectorSpace, 5, func(s []int) bool {
			b := []uint{}
			for _, idx := range s {
				b = append(b, fullVectorSpace[idx])
			}
			cl, err := codeloops.NewCL(codeloops.CLParams{Basis: b})
			if err != nil {
				log.Fatal(err)
			}
			if len(cl.VectorSpace()) != 1<<uint(len(s)) {
				return false
			}
			for _, v := range cl.VectorSpace() {
				if codeloops.BitWeight(v)%8 != 0 {
					return false
				}
			}
			return cl.IsAssoc()
		})
		if ok {
			log.Printf("Basis 1 found:")
			for _, idx := range set {
				b6 = append(b6, fullVectorSpace[idx])
				fmt.Printf("\t0x%x\n", fullVectorSpace[idx])
			}
			vs6 = append(vs6, codeloops.VectorSpace(b6)...)
			break
		}
	}

	for i := 0; i < len(fullVectorSpace); i++ {
		s := []int{i}
		ok, set := AddElems(s, fullVectorSpace, 4, func(s []int) bool {
			b := []uint{}
			for _, idx := range s {
				b = append(b, fullVectorSpace[idx])
			}
			cl, err := codeloops.NewCL(codeloops.CLParams{Basis: b})
			if err != nil {
				log.Fatal(err)
			}
			if len(cl.VectorSpace()) != 1<<uint(len(s)) {
				return false
			}
			for _, v := range cl.VectorSpace() {
				if codeloops.BitWeight(v)%8 != 0 {
					return false
				}
			}
			for _, v := range cl.VectorSpace() {
				for _, w := range vs6 {
					if v == w && v != 0 {
						return false
					}
				}
			}
			return cl.IsAssoc()
		})
		if ok {
			log.Printf("Basis 2 found:")
			for _, idx := range set {
				b5 = append(b5, fullVectorSpace[idx])
				fmt.Printf("\t0x%x\n", fullVectorSpace[idx])
			}
			vs5 = append(vs5, codeloops.VectorSpace(b5)...)
			break
		}
	}

	b11 := append(b6, b5...)
	vs11 := codeloops.VectorSpace(b11)
	if len(vs11) != 1<<11 {
		log.Fatalf("Not enough vectors in VS11: %d", len(vs11))
	}

	var lastVec uint
	for _, v := range fullVectorSpace {
		found := false
		for _, w := range vs11 {
			if v == w {
				found = true
				break
			}
		}
		if !found {
			log.Printf("Last vector found: 0x%x", v)
			lastVec = v
			break
		}
	}
	b12 := append(b11, lastVec)
	vs12 := codeloops.VectorSpace(b12)
	if len(vs12) != 1<<12 {
		log.Fatalf("Not enough vectors in VS11: %d", len(vs12))
	}

	log.Printf("Final Basis:")
	for _, v := range b12 {
		fmt.Printf("\t0x%x\n", v)
	}
}
