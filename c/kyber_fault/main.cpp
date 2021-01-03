//============================================================================
// Name        : kyber_fault.cpp
// Author      : 
// Version     :
// Copyright   : 
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <iomanip>
#include <cstdint>
#include <cstring>
#include <array>
#include <assert.h>
#include <random>
#include <climits>
#include <functional>
#include <fstream>
#include <sstream>
#include <unistd.h>

extern "C"{
#include "kyber/api.h"
#include "kyber/indcpa.h"
#include "kyber/randombytes.h"
#include "kyber/poly.h"
#include "kyber/symmetric.h"
#include "kyber/verify.h"

#if KYBER_K == 2 //fixed key is for kyber512
#include "fixed_keys.h" //fixed_keys is a c header, since we want to use the same file on the microcontroller
#endif
}

#include "basictests.h"
#include "maskeddec.h"

using namespace std;

//from rounding.cpp
extern void round_dist_v(int k=3); // @suppress("Unused function declaration")
extern void round_dist_u(int k=3); // @suppress("Unused function declaration")

// stores the used randombytes function
function<void(unsigned char*, size_t)> randombytes_func;

using random_bytes_engine = independent_bits_engine<mt19937, CHAR_BIT, unsigned char>;
random_bytes_engine rbe;
void randombytes_deterministic(unsigned char *buf,size_t buflen){
	for(auto i = 0u; i < buflen; i++){
		buf[i] = rbe();
	}
}

void openfile(ofstream& fs, string filename){
	fs.open(filename, ios::out | ios::trunc);
	if(!fs.is_open()){
		cerr << "Failed to open file " << filename << endl;
		exit(1);
	}
}

string hexstr(uint8_t* bytes, int n){
	stringstream ss;
	ss << hex;
	ss << setfill('0');
	for(auto i = 0; i < n; i++){
		ss << setw(2) << (unsigned)bytes[i];
	}
	return ss.str();
}

void polyvec_sub(polyvec *r, const polyvec *a, const polyvec *b)
{
  int i;
  for(i=0;i<KYBER_K;i++)
    poly_sub(&r->vec[i], &a->vec[i], &b->vec[i]);
}

//keypair generation + copying of the two noise vectors
void indcpa_keypair_raw(unsigned char *pk, unsigned char *sk, polyvec& key_s, polyvec& key_e)
{
  polyvec a[KYBER_K], e, pkpv, skpv;
  unsigned char buf[2*KYBER_SYMBYTES];
  unsigned char *publicseed = buf;
  unsigned char *noiseseed = buf+KYBER_SYMBYTES;
  int i;
  unsigned char nonce=0;

  randombytes_func(buf, KYBER_SYMBYTES);
  hash_g(buf, buf, KYBER_SYMBYTES);

  gen_a(a, publicseed);

  for(i=0;i<KYBER_K;i++)
    poly_getnoise(skpv.vec+i, noiseseed, nonce++);
  for(i=0;i<KYBER_K;i++)
    poly_getnoise(e.vec+i, noiseseed, nonce++);

  key_s = skpv;
  key_e = e;

  polyvec_ntt(&skpv);
  polyvec_ntt(&e);

  // matrix-vector multiplication
  for(i=0;i<KYBER_K;i++) {
    polyvec_pointwise_acc(&pkpv.vec[i], &a[i], &skpv);
    poly_frommont(&pkpv.vec[i]);
  }

  polyvec_add(&pkpv, &pkpv, &e);
  polyvec_reduce(&pkpv);

  pack_sk(sk, &skpv);
  pack_pk(pk, &pkpv, publicseed);
}

int crypto_kem_keypair_raw(unsigned char *pk, unsigned char *sk, polyvec& key_s, polyvec& key_e)
{
  size_t i;
  indcpa_keypair_raw(pk,sk, key_s, key_e);
  for(i=0;i<KYBER_INDCPA_PUBLICKEYBYTES;i++)
    sk[i+KYBER_INDCPA_SECRETKEYBYTES] = pk[i];
  hash_h(sk+KYBER_SECRETKEYBYTES-2*KYBER_SYMBYTES, pk, KYBER_PUBLICKEYBYTES);
  randombytes_func(sk+KYBER_SECRETKEYBYTES-KYBER_SYMBYTES, KYBER_SYMBYTES);        /* Value z for pseudo-random output on reject */
  return 0;
}

void indcpa_enc_raw(unsigned char *c, const unsigned char *m, const unsigned char *pk, const unsigned char *coins,
		polyvec& r, polyvec& e1, poly& e2,
		poly& v_rounded, poly& v_noround, poly& v_roundoff,
		polyvec& u_rounded, polyvec& u_noround, polyvec& u_roundoff)
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

  //copy noise polys before feeding them to the ntt
  r = sp;
  e1 = ep;
  e2 = epp;

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

  //copy ciphertext before rounding (in pack_ciphertext)
  v_noround = v;
  u_noround = bp;

  pack_ciphertext(c, &bp, &v);

  //we get the roundoff by unpacking again and checking the error
  unpack_ciphertext(&u_rounded, &v_rounded, c);
  //v_rounded = v_noround + v_roundoff --> v_roundoff = r_rounded - v_noround
  poly_sub(&v_roundoff, &v_rounded, &v_noround);
  for(auto &t: v_roundoff.coeffs){
	  if(t < -KYBER_Q/2){
		  t += KYBER_Q; //roundoff at zero
	  }
  }

  //u_rounded = u_noround + u_roundoff
  polyvec_sub(&u_roundoff, &u_rounded, &u_noround);
  for(auto& t1: u_roundoff.vec){
	  for(auto& t2: t1.coeffs){
		  if( t2 < -KYBER_Q/2)
			  t2 += KYBER_Q;
	  }
  }
}

int crypto_kem_enc_raw(unsigned char *ct, unsigned char *ss, const unsigned char *pk,
		polyvec& r, polyvec& e1, poly& e2,
		poly& v_rounded, poly& v_noround, poly& v_roundoff,
		polyvec& u_rounded, polyvec& u_noround, polyvec& u_roundoff,
		uint8_t* msg)
{
  unsigned char  kr[2*KYBER_SYMBYTES];                                     /* Will contain key, coins */
  unsigned char buf[2*KYBER_SYMBYTES];

  randombytes_func(buf, KYBER_SYMBYTES);
  hash_h(buf, buf, KYBER_SYMBYTES);                                        /* Don't release system RNG output */

  hash_h(buf+KYBER_SYMBYTES, pk, KYBER_PUBLICKEYBYTES);                    /* Multitarget countermeasure for coins + contributory KEM */
  hash_g(kr, buf, 2*KYBER_SYMBYTES);

  memcpy(msg, buf, KYBER_INDCPA_MSGBYTES);
  indcpa_enc_raw(ct, buf, pk, kr+KYBER_SYMBYTES, r, e1, e2, v_rounded, v_noround, v_roundoff, u_rounded, u_noround, u_roundoff);                              /* coins are in kr+KYBER_SYMBYTES */

  hash_h(kr+KYBER_SYMBYTES, ct, KYBER_CIPHERTEXTBYTES);                    /* overwrite coins in kr with H(c) */
  kdf(ss, kr, 2*KYBER_SYMBYTES);                                           /* hash concatenation of pre-k and H(c) to k */
  return 0;
}

void indcpa_dec_fault(unsigned char *m, const unsigned char *c, const unsigned char *sk, int faultcoeff)
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

//	poly_tomsg(m, &mp);
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
	    	if(8*i + j == faultcoeff)
	    		t = (((a->coeffs[8*i+j]) + KYBER_Q/2) / KYBER_Q) & 1;
	    	else
	    		t = (((a->coeffs[8*i+j] << 1) + KYBER_Q/2) / KYBER_Q) & 1;
	      msg[i] |= t << j;
	    }
	  }
	}
}

int crypto_kem_dec_fault(unsigned char *ss, const unsigned char *ct, const unsigned char *sk, int faultcoeff)
{
  size_t i;
  int fail;
  unsigned char __attribute__((aligned(32))) cmp[KYBER_CIPHERTEXTBYTES];
  unsigned char buf[2*KYBER_SYMBYTES];
  unsigned char kr[2*KYBER_SYMBYTES];                                      /* Will contain key, coins */
  const unsigned char *pk = sk+KYBER_INDCPA_SECRETKEYBYTES;

  indcpa_dec_fault(buf, ct, sk, faultcoeff);

  for(i=0;i<KYBER_SYMBYTES;i++)                                            /* Multitarget countermeasure for coins + contributory KEM */
    buf[KYBER_SYMBYTES+i] = sk[KYBER_SECRETKEYBYTES-2*KYBER_SYMBYTES+i];   /* Save hash by storing H(pk) in sk */
  hash_g(kr, buf, 2*KYBER_SYMBYTES);

  indcpa_enc(cmp, buf, pk, kr+KYBER_SYMBYTES);                             /* coins are in kr+KYBER_SYMBYTES */

  fail = verify(ct, cmp, KYBER_CIPHERTEXTBYTES);

  hash_h(kr+KYBER_SYMBYTES, ct, KYBER_CIPHERTEXTBYTES);                    /* overwrite coins in kr with H(c)  */

  cmov(kr, sk+KYBER_SECRETKEYBYTES-KYBER_SYMBYTES, KYBER_SYMBYTES, fail);  /* Overwrite pre-k with z on re-encryption failure */

  kdf(ss, kr, 2*KYBER_SYMBYTES);                                           /* hash concatenation of pre-k and H(c) to k */
  return 0;
}


// e*r - s*e_1
int16_t scalar_product(const polyvec& r, const polyvec& e_1, const polyvec& e, const polyvec& s, const int16_t coeff) {
	int16_t sum1 = 0;
	for (int16_t j = 0; j < KYBER_K; j++) {
		for (int16_t i = 0; i <= coeff; i++) {
			sum1 += r.vec[j].coeffs[coeff - i] * e.vec[j].coeffs[i];
		}
	}

	int16_t sum2 = 0;
	for (int16_t j = 0; j < KYBER_K; j++) {
		for (int16_t i = coeff + 1; i < KYBER_N; i++) {
			sum2 += r.vec[j].coeffs[KYBER_N + coeff - i] * e.vec[j].coeffs[i];
		}
	}

	int16_t sum3 = 0;
	for (int16_t j = 0; j < KYBER_K; j++) {
		for (int16_t i = 0; i <= coeff; i++) {
			sum3 += e_1.vec[j].coeffs[coeff - i] * s.vec[j].coeffs[i];
		}
	}

	int16_t sum4 = 0;
	for (int16_t j = 0; j < KYBER_K; j++) {
		for (int16_t i = coeff + 1; i < KYBER_N; i++) {
			sum4 += e_1.vec[j].coeffs[KYBER_N + coeff - i] * s.vec[j].coeffs[i];
		}
	}

	return sum1 - sum2 - sum3 + sum4;
}

string print_linear_coeffs(const polyvec& r, const polyvec& e_1, int16_t coeff) {
	stringstream ss;
	for (auto j = 0; j < KYBER_K; j++) {
		for (auto i = 0; i <= coeff; i++) {
			ss << setw(3) << r.vec[j].coeffs[coeff - i];
		}

		for (auto i = coeff + 1; i < KYBER_N; i++) {
			ss << setw(3) << -(r.vec[j].coeffs[KYBER_N + coeff - i]);
		}
	}
	for (auto j = 0; j < KYBER_K; j++) {
		for (auto i = 0; i <= coeff; i++) {
			ss << setw(3) << -(e_1.vec[j].coeffs[coeff - i]);
		}

		for (int16_t i = coeff + 1; i < KYBER_N; i++) {
			ss << setw(3) << e_1.vec[j].coeffs[KYBER_N + coeff - i];
		}
	}
	return ss.str();
}

void simple_hist(int numtests = 1000){

	array<uint64_t, KYBER_Q> hist_plain = {0};
	array<uint64_t, KYBER_Q> hist_err = {0};
	array<uint64_t, KYBER_Q> hist_v_roundoff = {0};
	array<uint64_t, KYBER_Q> hist_err_rcorrected = {0};
	array<uint64_t, KYBER_Q> hist_err_rcorrected2 = {0};


	uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES], sk[KYBER_INDCPA_SECRETKEYBYTES];
	uint8_t ct[KYBER_INDCPA_BYTES];
	uint8_t msg_a[KYBER_INDCPA_MSGBYTES], msg_b[KYBER_INDCPA_MSGBYTES];
	uint8_t coins[KYBER_SYMBYTES];

	for(int nt = 0; nt < numtests; nt++){
		randombytes_func(msg_a, KYBER_INDCPA_MSGBYTES);
		randombytes_func(coins, KYBER_SYMBYTES);

		polyvec key_s, key_e;
		indcpa_keypair_raw(pk, sk, key_s, key_e);

		//enc
		uint8_t* c = ct;
		uint8_t* m = msg_a;
		poly v_rounded, v_noround, v_roundoff;
		polyvec u_rounded, u_noround, u_roundoff;
		polyvec r, e1;
		poly e2;
		indcpa_enc_raw(c, m, pk, coins, r, e1, e2, v_rounded, v_noround, v_roundoff, u_rounded, u_noround, u_roundoff);
		for(auto &t: v_roundoff.coeffs){
			hist_v_roundoff[t + KYBER_Q/2]++;
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


				  for(i = 0; i < KYBER_N; i++){
					  hist_plain[a->coeffs[i]]++;
				  }

				  //re-encode known message and subtract
				  poly menc;
				  poly_frommsg(&menc, msg_a);
				  poly err;
				  poly_sub(&err, a, &menc);
				  for(auto& t: err.coeffs){
					  if(t > KYBER_Q/2)
						  t -= KYBER_Q;

					  hist_err[t+KYBER_Q/2]++;
				  }

				  //subtract roundoff
				  poly err_rcorrect;
				  poly_sub(&err_rcorrect, &err, &v_roundoff);
				  for(auto& t: err_rcorrect.coeffs){
					  hist_err_rcorrected[t+KYBER_Q/2]++;
				  }

				  //subtract e2
				  poly err_rcorrect2;
				  poly_sub(&err_rcorrect2, &err_rcorrect, &e2);
				  for(auto& t: err_rcorrect2.coeffs){
					  hist_err_rcorrected2[t+KYBER_Q/2]++;
				  }


				  for(i=0;i<KYBER_SYMBYTES;i++)
				  {
					msg[i] = 0;
					for(j=0;j<8;j++)
					{
					  t = (((a->coeffs[8*i+j] << 1) + KYBER_Q/2) / KYBER_Q) & 1;
					  msg[i] |= t << j;
					}
				  }

				  //sanity check: err_rcorrect2 should now be e*r - s*(e1 + u_roundoff)
				  {
					  polyvec_ntt(&key_s);
					  polyvec_ntt(&key_e);
					  polyvec_ntt(&r);

					  polyvec e1x;
					  polyvec_add(&e1x, &e1, &u_roundoff);
					  polyvec_ntt(&e1x);

					  poly acc1, acc2;
					  polyvec_pointwise_acc(&acc1, &key_e, &r);
					  polyvec_pointwise_acc(&acc2, &key_s, &e1x);
					  poly_invntt(&acc1);
					  poly_invntt(&acc2);

					  poly acc;
					  poly_sub(&acc, &acc1, &acc2);
					  poly_reduce(&acc);
					  poly_csubq(&acc);
					  for(auto& t: acc.coeffs){
						  if(t > KYBER_Q/2)
							  t -= KYBER_Q;
					  }

					  auto ok = memcmp(&acc, &err_rcorrect2, sizeof(poly)) == 0;
					  if(!ok){
						  cout << "FAIL error!!!" << endl;
					  }
				  }
			  }
		}

		auto eq = (memcmp(msg_a, msg_b, KYBER_INDCPA_MSGBYTES) == 0);
		if(!eq){
			cout << "FAIL msg!!" << endl;
			return;
		}
	}

	//print the final histogram
//	for(int i = 0; i < KYBER_Q; i++){
//		cout << setw(4) << i << ": " << hist[i] << endl;
//	}

//	for(int i = 0; i < KYBER_Q; i++){
//		cout << hist_plain[i] << endl;
//	}

//	for(int i = 0; i < KYBER_Q; i++){
//		cout << setw(5) << (i - KYBER_Q/2) << ": " << hist_err[i] << endl;
//	}

//	for(int i = 0; i < KYBER_Q; i++){
//		cout << setw(5) << (i - KYBER_Q/2) << " " << hist_v_roundoff[i] << endl;
//	}


//	for(int i = 0; i < KYBER_Q; i++){
//		cout << setw(5) << (i - KYBER_Q/2) << " " << hist_err_rcorrected[i] << endl;
//	}

	for(int i = 0; i < KYBER_Q; i++){
		cout << setw(5) << (i - KYBER_Q/2) << " " << hist_err_rcorrected2[i] << endl;
	}
}

void generate_system(int numfaults = 1000, int faultcoeff = 0, bool masked = false, int filter = -1){
	uint8_t pk[KYBER_INDCPA_PUBLICKEYBYTES], sk[KYBER_INDCPA_SECRETKEYBYTES];
	uint8_t ct[KYBER_INDCPA_BYTES];
	uint8_t msg_a[KYBER_INDCPA_MSGBYTES], msg_b[KYBER_INDCPA_MSGBYTES], msg_check[KYBER_INDCPA_MSGBYTES];
	uint8_t coins[KYBER_SYMBYTES];

	//generate a keypair and keep the two noise vectors, for sanity checking
	polyvec key_s, key_e;
	indcpa_keypair_raw(pk, sk, key_s, key_e);

	polyvec skpv1, skpv2;
	{ //unpack and share keys, do it in this block do prevent accidental use of skpv
		polyvec skpv;
		unpack_sk(&skpv, sk);
		mask_polyvec(skpv1, skpv2, skpv);
	}

	ofstream fs_key;
	openfile(fs_key, (masked ? "ineqs_es_masked.txt": "ineqs_es.txt"));

//	cout << "es: [";
	for(auto& p: key_e.vec){
		for(auto& px: p.coeffs){
			fs_key << setw(3) << px;
		}
	}
	for(auto& p: key_s.vec){
		for(auto& px: p.coeffs){
			fs_key << setw(3) << px;
		}
	}
//	cout << "]" << endl;
//	cout << endl;
	fs_key.close();

	ofstream fs_ineqs;
	openfile(fs_ineqs, (masked ? "ineqs_masked.txt" : "ineqs.txt"));

	//todo maybe multithreading to speed up the generation
	for(auto t = 0; t < numfaults; t++){
		randombytes_func(msg_a, KYBER_INDCPA_MSGBYTES);
		randombytes_func(coins, KYBER_SYMBYTES);

		// do the encryption
		poly v_rounded, v_noround, v_roundoff;
		polyvec u_rounded, u_noround, u_roundoff;
		polyvec r, e1;
		poly e2;
		indcpa_enc_raw(ct, msg_a, pk, coins, r, e1, e2, v_rounded, v_noround, v_roundoff, u_rounded, u_noround, u_roundoff);

		auto msgbit = (msg_a[faultcoeff/8] & (1 << (faultcoeff % 8))) != 0;

		//feed the ct to the faulty decryption
		if(!masked)
			indcpa_dec_fault(msg_b, ct, sk, faultcoeff);
		else
			indcpa_dec_fault_masked(msg_b, ct, &skpv1, &skpv2, faultcoeff);

		//test: is the whole message correct or not
		auto ineffective = (memcmp(msg_a, msg_b, KYBER_INDCPA_MSGBYTES) == 0);


		//ineffective: error term is positive ( >= 0)
		//effective: error term is negative ( < 0 )
		auto errorpositive = ineffective;

		//we then have the inequality
		// e*r - s (e_1 + du) (<, >=) -(e2 + dv)
		// e*l1 - s*l2 (<, >=) rhs

		auto rhs = -(e2.coeffs[faultcoeff] + v_roundoff.coeffs[faultcoeff]);
		if(masked)
			rhs -= (msgbit ? 2 : 1); // the masked decoder has some off-by-one or something, depending on the msgbit (-1 might be explained by > vs. >=, the second idk)

		if(filter >= 0 && abs(rhs) > filter){
			t--;
			continue;
		}

		auto l1 = r;
		polyvec l2;
		polyvec_add(&l2, &e1, &u_roundoff);

		//the masked decoder has a multiplication with -1 in there somewhere. Doesn't matter for decoding, but important to correct here
		if(masked){
			rhs = -rhs;
			polyvec_negate(&l1);
			polyvec_negate(&l2);
		}

		//print the inequalities
		fs_ineqs << print_linear_coeffs(l1, l2, faultcoeff);
		fs_ineqs << (errorpositive ? " >= " : " <  ") << setw(4) << rhs << endl;

		{
			//sanity check: compute the scalar product and check if inequality is correct
			auto sx = scalar_product(l1, l2, key_e, key_s, faultcoeff);

			auto check = (errorpositive ? sx >= rhs : sx < rhs);

			if(!check){
				cout << msgbit << endl;
				cout << "FAIL1!!!" << endl;
				return;
			}

			//also check if normal decryption still yields the correct result
			indcpa_dec(msg_check, ct, sk);
			check = (memcmp(msg_a, msg_check, KYBER_INDCPA_MSGBYTES) == 0);
			if(!check){
				cout << "FAIL2!!!" << endl;
				return;
			}
		}
	}
}

#if KYBER_K == 2
void generate_system_fixedkey(int numfaults = 1000, int faultcoeff = 0, int filter = -1){
	uint8_t *pk = fixed_pk;
	uint8_t *sk = fixed_sk;
	uint8_t ct[KYBER_CIPHERTEXTBYTES];
	uint8_t ss_a[KYBER_SSBYTES], ss_b[KYBER_SSBYTES];
	uint8_t msg[KYBER_INDCPA_BYTES];

	ofstream fs_ineqs;
	openfile(fs_ineqs, "ineqs_fixed.txt");

	//todo maybe multithreading to speed up the generation
	for(auto t = 0; t < numfaults; t++){
		// do the encryption
		poly v_rounded, v_noround, v_roundoff;
		polyvec u_rounded, u_noround, u_roundoff;
		polyvec r, e1;
		poly e2;

		crypto_kem_enc_raw(ct, ss_a, pk, r, e1, e2,  v_rounded, v_noround, v_roundoff, u_rounded, u_noround, u_roundoff, msg);

		//feed the ct to the faulty decryption
		crypto_kem_dec_fault(ss_b, ct, sk, faultcoeff);

		//test: is the whole message correct or not
		auto ineffective = (memcmp(ss_a, ss_b, KYBER_SSBYTES) == 0);


		//ineffective: error term is positive ( >= 0)
		//effective: error term is negative ( < 0 )
		auto errorpositive = ineffective;

		//we then have the inequality
		// e*r - s (e_1 + du) (<, >=) -(e2 + dv)
		// e*l1 - s*l2 (<, >=) rhs

		auto rhs = -(e2.coeffs[faultcoeff] + v_roundoff.coeffs[faultcoeff]);

		if(filter >= 0 && abs(rhs) > filter){
			t--;
			continue;
		}

		auto l1 = r;
		polyvec l2;
		polyvec_add(&l2, &e1, &u_roundoff);

		//print the inequalities
		fs_ineqs << print_linear_coeffs(l1, l2, faultcoeff);
		fs_ineqs << (errorpositive ? " >= " : " <  ") << setw(4) << rhs << endl;
	}
}

string generate_faultdata_single(int faultcoeff = 0, int filter = -1){
	uint8_t* pk = fixed_pk;
	uint8_t ct[KYBER_CIPHERTEXTBYTES];
	uint8_t ss[KYBER_SSBYTES];
	uint8_t msg[KYBER_INDCPA_MSGBYTES];

	stringstream ret;

	while(true){
		// do the encryption
		poly v_rounded, v_noround, v_roundoff;
		polyvec u_rounded, u_noround, u_roundoff;
		polyvec r, e1;
		poly e2;

		crypto_kem_enc_raw(ct, ss, pk, r, e1, e2, v_rounded, v_noround, v_roundoff, u_rounded, u_noround, u_roundoff, msg);

		//we then have the inequality
		// e*r - s (e_1 + du) (<, >=) -(e2 + dv)
		// e*l1 - s*l2 (<, >=) rhs

		auto rhs = -(e2.coeffs[faultcoeff] + v_roundoff.coeffs[faultcoeff]);

		if(filter >= 0 && abs(rhs) > filter){
			continue;
		}

		auto l1 = r;
		polyvec l2;
		polyvec_add(&l2, &e1, &u_roundoff);

		ret << hexstr(ct, KYBER_CIPHERTEXTBYTES) << ", ";
		ret << hexstr(ss, KYBER_SSBYTES) << ", ";
		ret << hexstr(msg, KYBER_INDCPA_MSGBYTES) << ", ";

		//print the inequalities
		ret << print_linear_coeffs(l1, l2, faultcoeff);
		ret <<  " ?? " << setw(4) << rhs;

		return ret.str();
	}
}

//generate data which will be used as input to the real fault injection
void generate_faultdata(int numfaults = 1000, int faultcoeff = 0, int filter = -1){
	ofstream fs;
	openfile(fs, "faultdata.txt");

	//todo maybe multithreading to speed up the generation
	for(auto t = 0; t < numfaults; t++){
		fs << generate_faultdata_single(faultcoeff, filter) << endl;
	}
	fs.close();
}

//generate a fixed key for the full kem (not indcpa-enc)
void gen_fixed_key(){
	// generate a keypair
	cout << "generate a keypair" << endl;
	uint8_t pk[KYBER_PUBLICKEYBYTES], sk[KYBER_SECRETKEYBYTES];
	polyvec key_s, key_e;
	crypto_kem_keypair_raw(pk, sk, key_s, key_e);

	// and write it to a file
	ofstream fs;
	openfile(fs, "fixed_keys.c");

	fs << "#include <stdint.h>" << endl << endl;

	fs << "const uint8_t fixed_pk[" << KYBER_PUBLICKEYBYTES <<  "] = {" << endl;

	fs << hex;
	fs << setfill('0');
	for(auto i = 0; i < KYBER_PUBLICKEYBYTES; i++){
		fs << "0x"<< setw(2) << (unsigned)pk[i]; //need to convert, otherwise treated as char and directly printed

		if(i < (KYBER_PUBLICKEYBYTES-1))
			fs << ", ";
		if((i+1)%10 == 0 || (i == KYBER_PUBLICKEYBYTES - 1))
			fs << endl;
	}
	fs << "};" << endl << endl;

	fs << dec;
	fs << "const uint8_t fixed_sk[" << KYBER_SECRETKEYBYTES <<  "] = {" << endl;
	fs << hex;
	for(auto i = 0; i < KYBER_SECRETKEYBYTES; i++){
		fs << "0x" << setw(2) << (unsigned)sk[i]; //need to convert, otherwise treated as char and directly printed

		if(i < (KYBER_SECRETKEYBYTES-1))
			fs << ", ";
		if((i+1)%10 == 0 || (i == KYBER_SECRETKEYBYTES - 1))
			fs << endl;
	}
	fs << "};" << endl;

	fs.close();

	openfile(fs, "fixed_keys.h");
	fs << dec;
	fs << "#ifndef FIXED_KEYS_H" << endl << "#define FIXED_KEYS_H" << endl << endl;
	fs << "#include <stdint.h>" << endl << endl;
	fs << "extern uint8_t fixed_pk[" << KYBER_PUBLICKEYBYTES << "];" << endl;
	fs << "extern uint8_t fixed_sk[" << KYBER_SECRETKEYBYTES << "];" << endl;

	fs << endl << endl << "#endif" << endl;

	fs.close();

	cout << "written keys to fixed_keys.c/h" << endl;

	// now write the key polynomials to a textfile, so that matlab/python understand them
	ofstream fs2; //no idea why, but with old fs the setw doesn't work anymore
	openfile(fs2, "fixed_es.txt");

	for(auto& p: key_e.vec){
		for(auto& px: p.coeffs){
			fs2 << setw(3) << px;
		}
	}
	for(auto& p: key_s.vec){
		for(auto& px: p.coeffs){
			fs2 << setw(3) << px;
		}
	}
	fs2 << endl;
	fs2.close();

	cout << "written sk-polys to es.txt" << endl;
}

//generate a ciphertext with the stored public key, print ciphertext and shared secret key
void gen_ct(){
	uint8_t ct[KYBER_CIPHERTEXTBYTES];
	uint8_t ss[KYBER_SSBYTES];
	crypto_kem_enc(ct, ss, fixed_pk);

	cout << hexstr(ct, KYBER_CIPHERTEXTBYTES) << endl;
	cout << hexstr(ss, KYBER_SSBYTES) << endl;
}

#endif


int main(int argc, char **argv) {

	//default parameters
	auto masked = false;
	auto seed = 0;
	rbe = random_bytes_engine(seed);
	randombytes_func = randombytes;
	auto numeqs = 10000;
	auto coeff = 1;
	auto filter = 10;

	int opt;

	// m: masked flag (default false)
	// s: seed (switches to deterministic mode)
	// n: number of equations
	// c: coefficient to attack
	// f: filter for ineqs, to only get those which have smaller than |f| on the rhs
	while((opt = getopt(argc, argv, "ms:n:c:f:")) != -1){
		switch(opt){
		case 'm':
			masked = true;
			break;
		case 's': //when a seed is specified, then we use deterministic generation
			seed = atoi(optarg);
			rbe = random_bytes_engine(seed);
			randombytes_func = randombytes_deterministic;
			break;
		case 'n':
			numeqs = atoi(optarg);
			break;
		case 'c':
			coeff = atoi(optarg);
			break;
		case 'f':
			filter = atoi(optarg);
			break;
		}
	}


	auto index = optind;
	if(argv[index] == NULL){ //standard
#if KYBER_K == 2
		cout << generate_faultdata_single(coeff, filter);
#else
		generate_system(numeqs, coeff, masked, filter);
#endif
		return 0;
	}

	//commands:
	// t: run testsuite
	// h: histogram
	// u: round_dist_u
	// v: round_dist_v
	// g: gen_fixed_key
	// c: gen_ct

	// s: generate system
	// a: generate a fault for the ARM cortex m

	auto cmd = argv[index][0];
	switch(cmd){
	case 't':
		basic_kem();
		basic_pke();
		basic_pke_raw();
		basic_kem_raw();
#if KYBER_K == 2
		test_fixed_key();
#endif
		test_masked();
		break;
	case 'h':
		simple_hist(numeqs);
		break;
	case 'u':
		round_dist_u(KYBER_K);
		break;
	case 'v':
		round_dist_v(KYBER_K);
		break;
#if KYBER_K == 2
	case 'g':
		gen_fixed_key();
		break;
	case 'c':
		gen_ct();
		break;
#endif

	case 's':
		generate_system(numeqs, coeff, masked, filter);
		break;
#if KYBER_K == 2
	case 'a':
		cout << generate_faultdata_single(coeff, filter);
		break;
#endif
	}

	return 0;
}
