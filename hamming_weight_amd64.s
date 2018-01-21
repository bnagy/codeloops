#include "textflag.h"

TEXT ·popCntq(SB),NOSPLIT,$0
	POPCNTQ    x+0(FP), AX
	MOVQ       AX, ret+8(FP)
	RET
