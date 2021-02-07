
#include "stdafx.h"

void bench_driver();

void ecc_driver();

int main(int argc, char* argv[]) {
  
  for (int i = 0; i < 20; ++i) 
  {
    uint512_t alice_private, alice_public;
    dh_gen_keypair(alice_private, alice_public);

    uint512_t bob_private, bob_public; 
    dh_gen_keypair(bob_private, bob_public);

    uint512_t alice_secret = dh_gen_secret(alice_private, bob_public);
    uint512_t bob_secret = dh_gen_secret(bob_private, alice_public);

    std::cout << "alice: private = 0x" << std::uppercase << std::hex << alice_private << std::endl 
      <<         "public         = 0x" << alice_public  << std::endl 
      <<         "secret         = 0x" << alice_secret  << std::endl << std::endl;

    std::cout << "bob  : private = 0x" << bob_private   << std::endl 
      <<         "public         = 0x" << bob_public    << std::endl 
      <<         "secret         = 0x" << bob_secret    << std::endl;

    std::cout << std::endl;

    if(alice_secret != bob_secret)
      std::cout << "alice_secret != bob_secret" << std::endl << std::endl;
  }

  if(isPrime(P))
    std::cout <<    "P = Prime      = 0x" << P;
  else std::cout << "P = Not Prime  = 0x" << P;

  std::cout << std::endl << std::endl;

  if (isPrime(G))
    std::cout << "G = Prime      = 0x" << G;
  else std::cout << "G = Not Prime  = 0x" << G;

  std::cout << std::endl << std::endl;
  
  bench_driver();

  std::cout << std::endl << std::endl;

  uint512_t bitcoin(0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f_cppui);
  bitcoin += 1;
  bitcoin /= 2;
  
  if (isPrime(bitcoin))
    std::cout << "bitcoin +1 /2 = Prime      = 0x" << bitcoin;
  else std::cout << "bitcoin +1 /2 = Not Prime  = 0x" << bitcoin;
  
  std::cout << std::endl;

  ecc_driver();

  return 0;
}

template <typename T>
T FiniteFieldElement<T>::Prime = T("0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f");

void ecc_driver()
{  
  typedef EllipticCurve<INT> ec_t;

  ec_t myEllipticCurve(-60, 17);

  // print out a little info and test some properties
  std::cout << std::endl << "The elliptic curve: \n" << std::dec << myEllipticCurve << std::endl;

  std::cout << std::endl << std::endl;

  auto G = myEllipticCurve.G1;

  mxws_t<INT> mxws_t;

  INT alicePrivKey = mxws_t(myEllipticCurve.Prime);

  INT bobPrivKey  = mxws_t(myEllipticCurve.Prime);

  auto start = std::chrono::steady_clock::now();
  
  auto alicePubKey = alicePrivKey * G;
  auto alicePubKeyCompressed = compress(alicePubKey); 

  auto end = std::chrono::steady_clock::now();

  auto total = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << "cost: ms per trial,      : " << total << std::endl << std::endl;
  
  auto bobPubKey   = bobPrivKey * G;
  auto bobPubKeyCompressed = compress(bobPubKey);
  
  std::cout << "alicePrivKey             : \n" << std::hex << alicePrivKey << std::endl;
  std::cout << "alicePubKey              : \n" << std::hex << alicePubKey << std::endl;
  std::cout << "Compressed alicePubKey   : \n" << alicePubKeyCompressed << std::endl;
  std::cout << "Decompressed alicePubKey : \n";
  std::cout << decompress(alicePubKeyCompressed, myEllipticCurve.Prime, myEllipticCurve) << std::endl << std::endl;
  
  std::cout << "bobPrivKey               : \n" << std::hex << bobPrivKey << std::endl;
  std::cout << "bobPubKey                : \n" << bobPubKey << std::endl;
  std::cout << "Compressed bobPubKey     : \n" << bobPubKeyCompressed << std::endl;
  std::cout << "Decompressed bobPubKey   : \n";
  std::cout << decompress(bobPubKeyCompressed, myEllipticCurve.Prime, myEllipticCurve) << std::endl;

  auto i = gcd(powm(bobPubKey.x().i_, 3, myEllipticCurve.Prime), myEllipticCurve.Prime);
  if ( i != 1)
    std::cout << "a and p are not coprime" << std::endl;
 
  std::cout << std::endl;

  auto aliceSharedKey = alicePrivKey * bobPubKey;
  auto bobSharedKey   = bobPrivKey   * alicePubKey;

  std::cout << "aliceSharedKey           : \n" << aliceSharedKey << std::endl;
  std::cout << "bobSharedKey             : \n" << bobSharedKey << std::endl << std::endl;

}

template <typename T>
void test(T n, T p) {
  struct Solution<T> sol = ts(n, p);
  printf("n = %llu\n", n);
  printf("p = %llu\n", p);
  if (sol.exists) {
    printf("root1 = %llu\n", sol.root1);
    printf("root2 = %llu\n", sol.root2);
  }
  else {
    printf("No solution exists\n");
  }
  printf("\n");
}

template <typename T>
void print_array(T* sbox)
{
for (int i = 0; i < sizeof(sbox); i++) {
  if (((i + 1) % 16) == 1)std::cout << std::endl;
  std::cout << std::setfill('0') << std::setw(2) << std::hex << int(sbox[i]) << " ";
}}