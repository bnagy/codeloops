package codeloops

const (
	m1q uint = 0x5555555555555555
	m2q      = 0x3333333333333333
	m4q      = 0x0f0f0f0f0f0f0f0f
	hq       = 0x0101010101010101
)

// lifted from https://github.com/steakknife/hamming
func BitWeight(x uint) uint {
	x -= (x >> 1) & m1q              // put count of each 2 bits into those 2 bits
	x = (x & m2q) + ((x >> 2) & m2q) // put count of each 4 bits into those 4 bits
	x = (x + (x >> 4)) & m4q         // put count of each 8 bits into those 8 bits
	return (x * hq) >> 56            // returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...
}

// based on https://filippo.io/callback-based-SetCombinationsWithoutReplacement-in-go/
func SetCombinationsWithReplacement(setSz, subsetSz uint, f func([]uint) error) {
	// For each combination of subsetSz elements out of setSz
	// call the function f passing a list of subsetSz integers in 0-setSz
	// with repetitions
	// If f returns an error then we bail early.
	indices := make([]uint, subsetSz)
	err := f(indices)
	if err != nil {
		return
	}
loop:
	for {
		var i int
		for i = int(subsetSz - 1); i >= 0; i-- {
			if indices[i] != setSz-1 {
				break
			}
		}
		if i < 0 {
			break
		}

		indices_i := indices[i]
		for k := i; k < int(subsetSz); k++ {
			indices[k] = indices_i + 1
		}
		err = f(indices)
		if err != nil {
			break loop
		}
	}
}

func SetCombinationsWithoutReplacement(setSz, subsetSz uint, f func([]uint)) {
	// For each combination of subsetSz elements out of setSz
	// call the function f passing a list of subsetSz integers in 0-setSz
	// without repetitions

	// TODO: switch to iterative algo
	s := make([]uint, subsetSz)
	last := subsetSz - 1
	var rc func(uint, uint)
	rc = func(i, next uint) {
		for j := next; j < setSz; j++ {
			s[i] = j
			if i == last {
				f(s)
			} else {
				rc(i+1, j+1)
			}
		}
		return
	}
	rc(0, 0)
}
