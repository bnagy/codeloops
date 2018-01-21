package codeloops

const (
	m1q uint = 0x5555555555555555
	m2q      = 0x3333333333333333
	m4q      = 0x0f0f0f0f0f0f0f0f
	hq       = 0x0101010101010101
)

// lifted from https://github.com/steakknife/hamming
func bitWeight(x uint) uint {
	x -= (x >> 1) & m1q              // put count of each 2 bits into those 2 bits
	x = (x & m2q) + ((x >> 2) & m2q) // put count of each 4 bits into those 4 bits
	x = (x + (x >> 4)) & m4q         // put count of each 8 bits into those 8 bits
	return (x * hq) >> 56            // returns left 8 bits of x + (x<<8) + (x<<16) + (x<<24) + ...
}
