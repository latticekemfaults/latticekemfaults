#ifndef INDCPA_H
#define INDCPA_H


#include "poly.h"
#include "polyvec.h"

// need to make this stuff public
void pack_pk(unsigned char *r, polyvec *pk, const unsigned char *seed);
void unpack_pk(polyvec *pk, unsigned char *seed, const unsigned char *packedpk);
void pack_sk(unsigned char *r, polyvec *sk);
void unpack_sk(polyvec *sk, const unsigned char *packedsk);
void pack_ciphertext(unsigned char *r, polyvec *b, poly *v);
void unpack_ciphertext(polyvec *b, poly *v, const unsigned char *c);
unsigned int rej_uniform(int16_t *r, unsigned int len, const unsigned char *buf, unsigned int buflen);
void gen_matrix(polyvec *a, const unsigned char *seed, int transposed);

#define gen_a(A,B)  gen_matrix(A,B,0)
#define gen_at(A,B) gen_matrix(A,B,1)


void indcpa_keypair(unsigned char *pk,
                    unsigned char *sk);

void indcpa_enc(unsigned char *c,
                const unsigned char *m,
                const unsigned char *pk,
                const unsigned char *coins);

void indcpa_dec(unsigned char *m,
                const unsigned char *c,
                const unsigned char *sk);


#endif
