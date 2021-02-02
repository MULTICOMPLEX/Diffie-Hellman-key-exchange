# A simple Diffie–Hellman key exchange algorithm

https://en.wikipedia.org/wiki/Diffie%E2%80%93Hellman_key_exchange

https://www.youtube.com/watch?v=bjWOG50PfdI

https://en.wikipedia.org/wiki/Primitive_root_modulo_n

https://en.wikipedia.org/wiki/Discrete_logarithm

Diffie–Hellman key exchange (DH) is a method of securely exchanging cryptographic keys over a public channel.

## How to use

1.Alice generates private&public key pair
```
    uint64_t alice_private, alice_public;
    dh_gen_keypair(alice_private, alice_public);
```
2.Bob generates private&public key pair
```
    uint64_t bob_private, bob_public;
    dh_gen_keypair(bob_private, bob_public);
```

3.Alice and Bob exchange their public key

4.Generate their own secret key by public key from the other

for Alice:
```
    uint64_t alice_secret = dh_gen_secret(alice_private, bob_public);
```
for Bob:
```
    uint64_t bob_secret = dh_gen_secret(bob_private, alice_public);
```
