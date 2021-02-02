
#include <iostream>
#include <vector>
#include "ecc.hpp"

INT Prime = 0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f_cppui;

typedef EllipticCurve<Prime> ec_t;

void ecc_driver()
{

	ec_t   myEllipticCurve(0, 7);


	// print out a little info and test some properties
	std::cout << "The elliptic curve: " << myEllipticCurve << "\n";

	std::cout << "\n\n";

	ec_t::Point P = myEllipticCurve[0];
	//std::cout << "some point P  = " << (312121212543534 * P) << "\n";

	INT privKey = 0x51897b64e85c3f714bba707e867914295a1377a7463a9dae8ea6a8b914246319_cppui;

	ec_t::Point pubKey = privKey * P;

	INT aa = pubKey.y_;

	std::cout << aa  << std::endl;

}
