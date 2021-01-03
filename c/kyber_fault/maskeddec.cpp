/*
 * maskeddec.cpp
 *
 *  Created on: Mar 16, 2020
 *      Author: ppessl
 */

#include "maskeddec.h"

#include <iostream>
#include <cstring>
#include <cstdint>
#include <functional>

extern "C"{
#include "kyber/indcpa.h"
#include "kyber/randombytes.h"
}

extern void indcpa_keypair_raw(unsigned char *pk, unsigned char *sk, polyvec& key_s, polyvec& key_e);

constexpr int16_t QM = (1 << QBITS) - 1;
constexpr int16_t QM2 = (1 << QBITS2) - 1;

using namespace std;
extern function<void(unsigned char*, size_t)> randombytes_func;

void poly_negate(poly* r){
	for(int i = 0u; i < KYBER_N; i++){
		r->coeffs[i] = -r->coeffs[i];
	}
//	poly_reduce(r);
}

void polyvec_negate(polyvec* r){
	for(int i = 0u; i < KYBER_K; i++){
		poly_negate(&(r->vec[i]));
	}
}

static int16_t rval(){
	int16_t ret;
	randombytes_func( (uint8_t*)&ret, 2);
	return (ret & QM);
}

static int16_t rval2(){
	int16_t ret;
	randombytes_func( (uint8_t*)&ret, 2);
	return (ret & QM2);
}

// generates a fully random polynomial
static void rand_poly(poly& r){
	for(auto& t: r.coeffs){
		do{
			t = rval();
		} while(t >= KYBER_Q);
	}
}

// generates a fully random polyvec
//static void rand_polyvec(polyvec& r){
//	for(auto& s: r.vec){
//		rand_poly(s);
//	}
//}

// generate a random sharing of s = s1+s2 mod q
void mask_poly(poly& s1, poly& s2, const poly& s){
	//first generate a random poly
	rand_poly(s1);

	//s2 = s - s1
	for(auto i = 0u; i < KYBER_N; i++){
		s2.coeffs[i] = (s.coeffs[i] - s1.coeffs[i] + KYBER_Q) % KYBER_Q;
	}
}


void mask_polyvec(polyvec& s1, polyvec& s2, const polyvec& s){
	for(auto i = 0u; i < KYBER_K; i++){
		mask_poly(s1.vec[i], s2.vec[i], s.vec[i]);
	}
}

void refresh_poly(poly& s1, poly& s2){
	poly r;
	rand_poly(r);

	poly_add(&s1, &s1, &r);
	poly_csubq(&s1);

	poly_sub(&s2, &s2, &r);
	for(auto& t: s2.coeffs){
		t = (t + KYBER_Q) % KYBER_Q;
	}
}

void refresh_polyvec(polyvec& s1, polyvec& s2){
	for(auto i = 0u; i < KYBER_K; i++){
		refresh_poly(s1.vec[i], s2.vec[i]);
	}
}

// arithmetic to boolean conversion
static void A2B(uint16_t& y1, uint16_t& y2, const uint16_t x1, const uint16_t x2){
	// TODO do a proper resharing, not via unmasking
	int16_t val = (x1 + x2) & QM2;

	y1 = rval2();
	y2 = y1 ^ val;
}

// from: Practical CCA2-Secure and Masked Ring-LWE Implementation
static void TransformPower2(uint16_t& y1, uint16_t& y2, const uint16_t x1, const uint16_t x2){
	y1 = rval2();
	y2 = x1 - y1;
	y2 = y2 + x2;
	uint16_t z1 = y1 - KYBER_Q;
	uint16_t z2 = 0;
	A2B(z1, z2, z1, y2);
	uint16_t k1 = MSB(z1) ^ 0x1;
	uint16_t k2 = MSB(z2);
	uint16_t k1p = rval2();
	uint16_t k1pp = k1 - k1p;
	uint16_t k2p = rval2();
	uint16_t k2pp = k2 - k2p;
	uint16_t r = rval2();
	y1 = (((((((r + y1) - k1*KYBER_Q) - k2*KYBER_Q) + 2*k1p*k2p*KYBER_Q) + 2*k1p*k2pp*KYBER_Q) + 2*k1pp*k2p*KYBER_Q) + 2*k1pp*k2pp*KYBER_Q) & QM2;
	y2 = (y2 - r) & QM2;
}

void masked_decode(uint16_t& m1, uint16_t& m2, int16_t c1, int16_t c2){
	// first step: -q/4 mod q
	c1 = (c1 - KYBER_Q/4 + KYBER_Q) % KYBER_Q;
	// second step: Transform
	uint16_t y1, y2;
	TransformPower2(y1, y2, c1, c2);
	// Third step: -q/2 mod 2^n
	y1 = (y1 - KYBER_Q/2) & QM2;
	// Fourth step: antoher a2b
	uint16_t y1p, y2p;
	A2B(y1p, y2p, y1, y2);
	// Fifth step: extract MSB, done
	m1 = MSB(y1p);
	m2 = MSB(y2p);
}

void masked_decode_fault(uint16_t& m1, uint16_t& m2, int16_t c1, int16_t c2){
	// first step: -q/4 mod q
//	c1 = (c1 - KYBER_Q/4 + KYBER_Q) % KYBER_Q; // skip the subtraction by q/4
	// second step: Transform
	uint16_t y1, y2;
	TransformPower2(y1, y2, c1, c2);
	// Third step: -q/2 mod 2^n
	y1 = (y1 - KYBER_Q/2) & QM2;
	// Fourth step: antoher a2b
	uint16_t y1p, y2p;
	A2B(y1p, y2p, y1, y2);
	// Fifth step: extract MSB, done
	m1 = MSB(y1p);
	m2 = MSB(y2p);
}

void test_masked(){
	uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES], sk[KYBER_INDCPA_SECRETKEYBYTES];
	uint8_t ct[KYBER_INDCPA_BYTES];
	uint8_t msg_a[KYBER_INDCPA_MSGBYTES], msg_b[KYBER_INDCPA_MSGBYTES];
	uint8_t coins[KYBER_SYMBYTES];

	randombytes_func(msg_a, KYBER_INDCPA_MSGBYTES);
	randombytes_func(coins, KYBER_SYMBYTES);

	polyvec key_s, key_e;
	indcpa_keypair_raw(pk, sk, key_s, key_e); //need to use the raw version to also allow deterministic behavior
	polyvec skpv1, skpv2;
	{ //unpack and share keys, do it in this block do prevent accidental use of skpv
		polyvec skpv;
		unpack_sk(&skpv, sk);
		mask_polyvec(skpv1, skpv2, skpv);
	}

	indcpa_enc(ct, msg_a, pk, coins);

//	{
//		polyvec bp;
//		poly v, mp1, mp2;
//
//		unpack_ciphertext(&bp, &v, ct);
//		polyvec_ntt(&bp);
//
//		//process first share
//		polyvec_pointwise_acc(&mp1, &skpv1, &bp);
//		poly_invntt(&mp1);
//
//		poly_sub(&mp1, &v, &mp1);
//		poly_reduce(&mp1);
//
//		//process second share
//		polyvec_pointwise_acc(&mp2, &skpv2, &bp);
//		poly_invntt(&mp2);
//		poly_negate(&mp2); //no subtraction here, but need to negate (v - u*s)
//
//
//	//	poly_tomsg(m, &mp);
//		{
//		  uint16_t t;
//		  int i,j;
////		  poly *a = &mp;
//		  uint8_t* msg = msg_b;
//
//		  for(i=0;i<KYBER_SYMBYTES;i++)
//		  {
//		    msg[i] = 0;
//		    for(j=0;j<8;j++)
//		    {
////		    	t = (((mp.coeffs[8*i+j] << 1) + KYBER_Q/2) / KYBER_Q) & 1;
//		    	uint16_t t0, t1;
//		    	masked_decode(t0, t1, mp1.coeffs[8*i + j], mp2.coeffs[8*i + j]);
//		    	t = t0^t1; //unmask here
//		      msg[i] |= t << j;
//		    }
//		  }
//		}
//	}

	indcpa_dec_fault_masked(msg_b, ct, &skpv1, &skpv2, -1);

	auto eq = (memcmp(msg_a, msg_b, KYBER_INDCPA_MSGBYTES) == 0);
	cout << "Masked PKE: " <<  (eq ? "Success" : "Fail") << endl;
}


void indcpa_dec_fault_masked(unsigned char *m, const unsigned char *ct, polyvec* skpv1, polyvec* skpv2, int faultcoeff){
	polyvec bp;
	poly v, mp1, mp2;

	unpack_ciphertext(&bp, &v, ct);
	polyvec_ntt(&bp);

	//process first share
	polyvec_pointwise_acc(&mp1, skpv1, &bp);
	poly_invntt(&mp1);

	poly_sub(&mp1, &v, &mp1);
	poly_reduce(&mp1);

	//process second share
	polyvec_pointwise_acc(&mp2, skpv2, &bp);
	poly_invntt(&mp2);
	poly_negate(&mp2); //no subtraction here, but need to negate (v - u*s)
	poly_reduce(&mp2);


	//	poly_tomsg(m, &mp);
	{
		uint16_t t;
		int i,j;

		for(i=0;i<KYBER_SYMBYTES;i++)
		{
			m[i] = 0;
			for(j=0;j<8;j++)
			{
				//		    	t = (((mp.coeffs[8*i+j] << 1) + KYBER_Q/2) / KYBER_Q) & 1;
				uint16_t t0, t1;
				if(8*i + j == faultcoeff)
					masked_decode_fault(t0, t1, mp1.coeffs[8*i + j], mp2.coeffs[8*i + j]);
				else
					masked_decode(t0, t1, mp1.coeffs[8*i + j], mp2.coeffs[8*i + j]);
				t = t0^t1; //unmask here
				m[i] |= t << j;
			}
		}
	}
}
