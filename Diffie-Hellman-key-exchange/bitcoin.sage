

E = EllipticCurve(GF(0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f, 'x'), [0, 7])
P = E.gen(0);P

B = E(3156566521962345705282722233281194109376473116309193527805390219356486815894,
     50621615428923838952367799849558827718095819354153345129785555578381658988628)
	 

B = F(3156566521962345705282722233281194109376473116309193527805390219356486815894,
     50621615428923838952367799849558827718095819354153345129785555578381658988628)
P*B

hex(P.order())

#
def rdp(nbits=256):
    while True:
        p = random_prime(2^nbits-1, false, 2^(nbits-1))
        if ZZ((p+1)/2).is_prime():
            return p

p = hex(rdp())

E = EllipticCurve(GF(p, 'x'), [0, 7])

P = E.gen(0);P
hex(P.order())
#

