package codeloops

// These return int because that's the natural type for comparison literals.

func bitWeight(x uint) (weight int) {
	// One loop iteration per set bit
	// https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetKernighan
	for weight = 0; x != 0; weight++ {
		x &= x - 1 // clear the least significant bit set
	}
	return
}

func parity(x uint) int {
	return bitWeight(x) & 1
}
