/*
 * basictests.cpp
 *
 *  Created on: Mar 17, 2020
 *      Author: ppessl
 */


#include <iostream>
#include <cstdint>
#include <cstring>

extern "C"{
#include "kyber/api.h"
#include "kyber/indcpa.h"
#include "kyber/randombytes.h"
#include "kyber/poly.h"

#include "fixed_keys.h" //fixed_keys is a c header, since we want to use the same file on the microcontroller
}

using namespace std;

void basic_kem(){
	uint8_t pk[KYBER_PUBLICKEYBYTES], sk[KYBER_SECRETKEYBYTES];
	uint8_t ct[KYBER_CIPHERTEXTBYTES];
	uint8_t ss_a[KYBER_SSBYTES], ss_b[KYBER_SSBYTES];

	crypto_kem_keypair(pk, sk);
	crypto_kem_enc(ct, ss_a, pk);
	crypto_kem_dec(ss_b, ct, sk);

	auto eq = (memcmp(ss_a, ss_b, KYBER_SSBYTES) == 0);

	cout << "Basic KEM: " <<  (eq ? "Success" : "Fail") << endl;
}

void basic_pke(){
	uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES], sk[KYBER_INDCPA_SECRETKEYBYTES];
	uint8_t ct[KYBER_INDCPA_BYTES];
	uint8_t msg_a[KYBER_INDCPA_MSGBYTES], msg_b[KYBER_INDCPA_MSGBYTES];
	uint8_t coins[KYBER_SYMBYTES];

	randombytes(msg_a, KYBER_INDCPA_MSGBYTES);
	randombytes(coins, KYBER_SYMBYTES);

	indcpa_keypair(pk, sk);
	indcpa_enc(ct, msg_a, pk, coins);
	indcpa_dec(msg_b, ct, sk);

	auto eq = (memcmp(msg_a, msg_b, KYBER_INDCPA_MSGBYTES) == 0);
	cout << "Basic PKE: " <<  (eq ? "Success" : "Fail") << endl;
}

extern void indcpa_keypair_raw(unsigned char *pk, unsigned char *sk, polyvec& key_s, polyvec& key_e);
void basic_pke_raw(){
	uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES], sk[KYBER_INDCPA_SECRETKEYBYTES];
	uint8_t ct[KYBER_INDCPA_BYTES];
	uint8_t msg_a[KYBER_INDCPA_MSGBYTES], msg_b[KYBER_INDCPA_MSGBYTES];
	uint8_t coins[KYBER_SYMBYTES];

	randombytes(msg_a, KYBER_INDCPA_MSGBYTES);
	randombytes(coins, KYBER_SYMBYTES);

	polyvec key_e, key_s;
	indcpa_keypair_raw(pk, sk, key_s, key_e);

	//enc
	uint8_t* c = ct;
	uint8_t* m = msg_a;
	{
		  polyvec sp, pkpv, ep, at[KYBER_K], bp;
		  poly v, k, epp;
		  unsigned char seed[KYBER_SYMBYTES];
		  int i;
		  unsigned char nonce=0;

		  unpack_pk(&pkpv, seed, pk);
		  poly_frommsg(&k, m);
		  gen_at(at, seed);

		  for(i=0;i<KYBER_K;i++)
		    poly_getnoise(sp.vec+i, coins, nonce++);
		  for(i=0;i<KYBER_K;i++)
		    poly_getnoise(ep.vec+i, coins, nonce++);
		  poly_getnoise(&epp, coins, nonce++);

		  polyvec_ntt(&sp);

		  // matrix-vector multiplication
		  for(i=0;i<KYBER_K;i++)
		    polyvec_pointwise_acc(&bp.vec[i], &at[i], &sp);

		  polyvec_pointwise_acc(&v, &pkpv, &sp);

		  polyvec_invntt(&bp);
		  poly_invntt(&v);

		  polyvec_add(&bp, &bp, &ep);
		  poly_add(&v, &v, &epp);
		  poly_add(&v, &v, &k);
		  polyvec_reduce(&bp);
		  poly_reduce(&v);

		  pack_ciphertext(c, &bp, &v);
	}

	//dec
	m = msg_b;
	{
		  polyvec bp, skpv;
		  poly v, mp;

		  unpack_ciphertext(&bp, &v, c);
		  unpack_sk(&skpv, sk);

		  polyvec_ntt(&bp);
		  polyvec_pointwise_acc(&mp, &skpv, &bp);
		  poly_invntt(&mp);

		  poly_sub(&mp, &v, &mp);
		  poly_reduce(&mp);

//		  poly_tomsg(m, &mp);
		  {
			  uint16_t t;
			  int i,j;
			  poly *a = &mp;
			  uint8_t* msg = m;

			  poly_csubq(a);


			  for(i=0;i<KYBER_SYMBYTES;i++)
			  {
			    msg[i] = 0;
			    for(j=0;j<8;j++)
			    {
			      t = (((a->coeffs[8*i+j] << 1) + KYBER_Q/2) / KYBER_Q) & 1;
			      msg[i] |= t << j;
			    }
			  }
		  }
	}

	auto eq = (memcmp(msg_a, msg_b, KYBER_INDCPA_MSGBYTES) == 0);
	cout << "Basic PKE raw: " <<  (eq ? "Success" : "Fail") << endl;
}

extern int crypto_kem_enc_raw(unsigned char *ct, unsigned char *ss, const unsigned char *pk,
		polyvec& r, polyvec& e1, poly& e2,
		poly& v_rounded, poly& v_noround, poly& v_roundoff,
		polyvec& u_rounded, polyvec& u_noround, polyvec& u_roundoff, uint8_t* msg);
void basic_kem_raw(){
	uint8_t pk[KYBER_PUBLICKEYBYTES], sk[KYBER_SECRETKEYBYTES];
	uint8_t ct[KYBER_CIPHERTEXTBYTES];
	uint8_t ss_a[KYBER_SSBYTES], ss_b[KYBER_SSBYTES];
	uint8_t msg[KYBER_INDCPA_MSGBYTES];


	polyvec r, e1;
	poly e2;
	poly v_rounded, v_noround, v_roundoff;
	polyvec u_rounded, u_noround, u_roundoff;

	{
		crypto_kem_keypair(pk, sk);
		crypto_kem_enc_raw(ct, ss_a, pk, r, e1, e2, v_rounded, v_noround, v_roundoff, u_rounded, u_noround, u_roundoff, msg);
		crypto_kem_dec(ss_b, ct, sk);

		auto eq = (memcmp(ss_a, ss_b, KYBER_SSBYTES) == 0);

		cout << "Basic KEM raw 1: " <<  (eq ? "Success" : "Fail") << endl;
	}

	{
		crypto_kem_enc_raw(ct, ss_a, fixed_pk, r, e1, e2, v_rounded, v_noround, v_roundoff, u_rounded, u_noround, u_roundoff, msg);
		crypto_kem_dec(ss_b, ct, fixed_sk);

		auto eq = (memcmp(ss_a, ss_b, KYBER_SSBYTES) == 0);

		cout << "Basic KEM raw 2: " <<  (eq ? "Success" : "Fail") << endl;
	}
}

#if KYBER_K == 2
void test_fixed_key(){
	uint8_t ct[KYBER_CIPHERTEXTBYTES];
	uint8_t ss_a[KYBER_SSBYTES], ss_b[KYBER_SSBYTES];
	crypto_kem_enc(ct, ss_a, fixed_pk);
	crypto_kem_dec(ss_b, ct, fixed_sk);

	auto eq = (memcmp(ss_a, ss_b, KYBER_SSBYTES) == 0);

	cout << "Stored-Key KEM: " <<  (eq ? "Success" : "Fail") << endl;
}
#endif
