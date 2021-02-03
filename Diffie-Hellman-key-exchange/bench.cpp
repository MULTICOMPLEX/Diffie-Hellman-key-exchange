
#include "dh.hpp"

std::pair<uint512_t, uint512_t> bench(unsigned int times);

void bench_driver() 
{
    using namespace std;
    const unsigned int times = 100;
    auto start = std::chrono::steady_clock::now();
    bench(times);
    auto end = std::chrono::steady_clock::now();
    auto total = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout << "cost: us per trial,          " << (double)total / (double)times << std::endl;

    std::cout << std::endl;

    start = std::chrono::steady_clock::now();
    rand_t<uint512_t>();
    end = std::chrono::steady_clock::now();
    total = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    std::cout << "cost: ns rand_t(),       " << std::dec << total << std::endl;
    
    start = std::chrono::steady_clock::now();
    rand_t_prime<uint512_t>();
    end = std::chrono::steady_clock::now();
    total = chrono::duration_cast<chrono::microseconds>(end - start).count();
    std::cout << "cost: us rand_t_prime(), " << total << std::endl;
}

std::pair<uint512_t, uint512_t> bench(unsigned int times)
{
  std::pair<uint512_t, uint512_t> p;

  uint512_t alice_private, alice_public;
  uint512_t bob_private, bob_public;

  uint512_t alice_secret = 0;
  uint512_t bob_secret = 0;

  for (unsigned int i = 0; i < times; ++i)
  {
    dh_gen_keypair(alice_private, alice_public);

    dh_gen_keypair(bob_private, bob_public);

    alice_secret = dh_gen_secret(alice_private, bob_public);
    bob_secret = dh_gen_secret(bob_private, alice_public);
  }

  p.first = alice_secret;
  p.second = bob_secret;

  return p;
}