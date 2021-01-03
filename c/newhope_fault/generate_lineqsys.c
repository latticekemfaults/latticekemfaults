#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <stdbool.h>

#include "cpapke.h"
#include "fips202.h"
#include "poly.h"
#include "rng.h"
#include "api.h"

#define ES_FILE "ineqs_es.txt"
#define RE1_FILE "ineqs.txt"

int16_t center(uint16_t x){
	int16_t v = x % NEWHOPE_Q;
	if(v > NEWHOPE_Q/2)
		v -= NEWHOPE_Q;
	return v;
}

void print_coeff(FILE *fp, uint16_t coeff, bool flip)
{
  fprintf(fp, " ");
  int16_t c = center(coeff);
  if(flip)
	  c = -c;
  fprintf(fp, "%d", c);
}



// <copied-from file="poly.c" function="coeff_freeze">
static uint16_t coeff_freeze(uint16_t x)
{
  uint16_t m,r;
  int16_t c;
  r = x % NEWHOPE_Q;

  m = r - NEWHOPE_Q;
  c = m;
  c >>= 15;
  r = m ^ ((r^m)&c);

  return r;
}
// </copied-from>

// <copied-from file="poly.c" function="flipabs">
static uint16_t flipabs(uint16_t x)
{
  int16_t r,m;
  r = coeff_freeze(x);

  r = r - NEWHOPE_Q/2;
  m = r >> 15;
  return (r + m) ^ m;
}
// </copied-from>

// <inspired-by file="poly.c" function="flipabs">
static uint16_t flipabs_faulted(uint16_t x)
{
  int16_t r,m;
  r = coeff_freeze(x);

  r = r - NEWHOPE_Q/2;
  m = r >> 15;
  // <fault>
  return (r + m);
  // </fault>
}
// </inspired-by>

// <copied-from file="cpapke.c" function="decode_c">
static void decode_c(poly *b, poly *v, const unsigned char *r)
{
  poly_frombytes(b, r);
  poly_decompress(v, r+NEWHOPE_POLYBYTES);
}
// </copied-from>

// <copied-from file="cpapke.c" function="decode_pk">
static void decode_pk(poly *pk, unsigned char *seed, const unsigned char *r)
{
  int i;
  poly_frombytes(pk, r);
  for(i=0;i<NEWHOPE_SYMBYTES;i++)
    seed[i] = r[NEWHOPE_POLYBYTES+i];
}
// </copied-from>

// <copied-from file="cpapke.c" function="gen_a">
static void gen_a(poly *a, const unsigned char *seed)
{
  poly_uniform(a,seed);
}
// </copied-from>

// <copied-from file="cpapke.c" function="encode_c">
static void encode_c(unsigned char *r, const poly *b, const poly *v)
{
  poly_tobytes(r,b);
  poly_compress(r+NEWHOPE_POLYBYTES,v);
}
// </copied-from>

// <copied-from file="cpapke.c" function="encode_pk">
static void encode_pk(unsigned char *r, const poly *pk, const unsigned char *seed)
{
  int i;
  poly_tobytes(r, pk);
  for(i=0;i<NEWHOPE_SYMBYTES;i++)
    r[NEWHOPE_POLYBYTES+i] = seed[i];
}
// </copied-from>

// Compute the scalar product “r⋅e - e_1⋅s” for the given target coefficient
int16_t
scalar_product(const poly *r, const poly *e_1, const poly *e, const poly *s, const int16_t coeff)
{
  int16_t sum1 = 0;
  for (int16_t i = 0; i < coeff + 1; i++)
    sum1 += center(r->coeffs[coeff - i]) * center(e->coeffs[i]);

  int16_t sum2 = 0;
  for (int16_t i = coeff + 1; i < NEWHOPE_N; i++)
    sum2 += center(r->coeffs[NEWHOPE_N + coeff - i]) * center(e->coeffs[i]);

  int16_t sum3 = 0;
  for (int16_t i = 0; i < coeff + 1; i++)
    sum3 += center(e_1->coeffs[coeff - i]) * center(s->coeffs[i]);

  int16_t sum4 = 0;
  for (int16_t i = coeff + 1; i < NEWHOPE_N; i++)
    sum4 += center(e_1->coeffs[NEWHOPE_N + coeff - i]) * center(s->coeffs[i]);

  return sum1 - sum2 - sum3 + sum4;
}

// Compute the scalar product “r⋅e - e_1⋅s” for the given target coefficient
void
print_linear_coeffs(FILE *fd, const poly *r, const poly *e_1, const int16_t coeff)
{
  for (int16_t i = 0; i < coeff + 1; i++)
    print_coeff(fd, r->coeffs[coeff - i], false);

  for (int16_t i = coeff + 1; i < NEWHOPE_N; i++)
    print_coeff(fd, r->coeffs[NEWHOPE_N + coeff - i], true);

  for (int16_t i = 0; i < coeff + 1; i++)
    print_coeff(fd, e_1->coeffs[coeff - i], true);

  for (int16_t i = coeff + 1; i < NEWHOPE_N; i++)
    print_coeff(fd, e_1->coeffs[NEWHOPE_N + coeff - i], false);
}

void
fprintBstr(FILE *fp, char *S, unsigned char *A, unsigned long long L)
{
  unsigned long long  i;
  fprintf(fp, "%s", S);

  for ( i=0; i<L; i++ )
    fprintf(fp, "%02X", A[i]);
  if ( L == 0 )
    fprintf(fp, "00");

  fprintf(fp, "\n");
}

void
fprintPoly(FILE *fp, char *S, poly *p)
{
  unsigned long long  i;
  fprintf(fp, "%s", S);

  fprintf(fp, "[");
  for ( i=0; i<NEWHOPE_N; i++ )
    if ( i == NEWHOPE_N - 1 )
      fprintf(fp, "%d", center(p->coeffs[i] % NEWHOPE_Q));
    else
      fprintf(fp, "%d ", center(p->coeffs[i] % NEWHOPE_Q));
  fprintf(fp, "]");

  fprintf(fp, "\n");
}

int generate_system(int numfaults, int target_coefficient, int filter){
  unsigned char msg[CRYPTO_BYTES];

  unsigned char k_coins_d[3*NEWHOPE_SYMBYTES];
  randombytes(k_coins_d, 3*NEWHOPE_SYMBYTES);

  FILE *ineqs_es = fopen(ES_FILE, "w");
  if (ineqs_es == NULL) {
	printf("file '%s' could not be created.", ES_FILE);
	return 2;
  }

  FILE *ineqs_re1 = fopen(RE1_FILE, "w");
  if (ineqs_re1 == NULL) {
	printf("file '%s' could not be created.", RE1_FILE);
	return 2;
  }

  // ================= ===================
  // LPR notation      NewHope ref impl
  // ----------------- -------------------
  // e                 invntt(ehat)
  // r                 invntt(sprime)
  // e_1               invntt(eprime)
  // s                 invntt(shat)
  // e_2               eprimeprime
  // ================= ===================
  poly e, r, e_1, s, e_2;

  // Generate the public/private keypair
  unsigned char pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];

  // <inspired-by file="cpapke.c" function="cpapke_keypair">
  // void cpapke_keypair(unsigned char *pk, unsigned char *sk)
  {
	poly ahat, ehat, ahat_shat, bhat, shat;
	unsigned char z[2*NEWHOPE_SYMBYTES];
	unsigned char *publicseed = z;
	unsigned char *noiseseed = z+NEWHOPE_SYMBYTES;

	randombytes(z, NEWHOPE_SYMBYTES);
	shake256(z, 2*NEWHOPE_SYMBYTES, z, NEWHOPE_SYMBYTES);

	gen_a(&ahat, publicseed);

	poly_sample(&shat, noiseseed, 0);
	poly_ntt(&shat);

	poly_sample(&ehat, noiseseed, 1);
	poly_ntt(&ehat);

	poly_mul_pointwise(&ahat_shat, &shat, &ahat);
	poly_add(&bhat, &ehat, &ahat_shat);

	poly_tobytes(sk, &shat);
	encode_pk(pk, &bhat, publicseed);

	// <save-data>
	e = ehat;
	s = shat;
	// </save-data>
  }
  // </inspired-by>

  // write ES_FILE
  poly_invntt(&e);
  poly_invntt(&s);

  for (int i = 0; i < NEWHOPE_N; i++)
	print_coeff(ineqs_es, e.coeffs[i], 0);
  for (int i = 0; i < NEWHOPE_N; i++)
	print_coeff(ineqs_es, s.coeffs[i], 0);
  fprintf(ineqs_es, "\n");

  // Encrypt message msg
  for (int t = 0; t < numfaults; t++)
  {
	unsigned char ct[CRYPTO_CIPHERTEXTBYTES];

	poly compression_error;

	// for each fault we take a different message, s', e' and e''
	randombytes(k_coins_d, 3*NEWHOPE_SYMBYTES);
	memset(msg, 0, NEWHOPE_SYMBYTES);

	// cpapke_enc(ct, msg, pk, k_coins_d+NEWHOPE_SYMBYTES);
	{
	  poly enc_sprime, enc_eprime, enc_vprime, enc_ahat, enc_bhat, enc_eprimeprime, enc_uhat, enc_v;
	  unsigned char publicseed[NEWHOPE_SYMBYTES];

	  poly_frommsg(&enc_v, msg);

	  decode_pk(&enc_bhat, publicseed, pk);
	  gen_a(&enc_ahat, publicseed);

	  poly_sample(&enc_sprime, k_coins_d+NEWHOPE_SYMBYTES, 0);
	  poly_sample(&enc_eprime, k_coins_d+NEWHOPE_SYMBYTES, 1);
	  poly_sample(&enc_eprimeprime, k_coins_d+NEWHOPE_SYMBYTES, 2);

	  poly_ntt(&enc_sprime);
	  poly_ntt(&enc_eprime);

	  poly_mul_pointwise(&enc_uhat, &enc_ahat, &enc_sprime);
	  poly_add(&enc_uhat, &enc_uhat, &enc_eprime);

	  poly_mul_pointwise(&enc_vprime, &enc_bhat, &enc_sprime);
	  poly_invntt(&enc_vprime);

	  poly_add(&enc_vprime, &enc_vprime, &enc_eprimeprime);
	  poly_add(&enc_vprime, &enc_vprime, &enc_v); // add message

	  // <save-data>
	  compression_error = enc_vprime;

	  r = enc_sprime;
	  e_1 = enc_eprime;
	  e_2 = enc_eprimeprime;
	  // </save-data>

	  encode_c(ct, &enc_uhat, &enc_vprime);
	}

	poly_invntt(&r);
	poly_invntt(&e_1);

	// <evaluate object="compression error">
	poly decompressed;
	poly_decompress(&decompressed, ct+NEWHOPE_POLYBYTES);
	for (int i = 0; i < NEWHOPE_N; i++) {
	  // REMINDER compression_error contains coefficients of enc_vprime
	  /// ⇒ error := - vprime + decompressed

	  compression_error.coeffs[i] = (-compression_error.coeffs[i] + decompressed.coeffs[i] + NEWHOPE_Q) % NEWHOPE_Q;
	}
	// sanity check: compression creates 8 uniformly-sized buckets, each of size 12289/8,
	//  thus compression error cannot exceed 1537
//    for (int i = 0; i < NEWHOPE_N; i++) {
//      int16_t val = center(compression_error.coeffs[i]);
//      if (abs(val) > NEWHOPE_Q/8 && !(decompressed.coeffs[i] == 0 && val >= 10753 && val < NEWHOPE_Q)) {
//        printf("FAIL compression error abs(%d) exceeds %u\n", val, NEWHOPE_Q/8);
//        return 5;
//      }
//    }
	// </evaluate>

	poly_add(&e_2, &e_2, &compression_error);
	int16_t rhs = -center(e_2.coeffs[target_coefficient] % NEWHOPE_Q) + 1;

	if(filter >= 0 && abs(rhs) > filter){
		t--;
		continue;
	}

	// Decrypt ciphertext
	unsigned char out[CRYPTO_BYTES];

	//cpapke_dec(out, ct, sk);
	{
	  // <inspired-by file="cpapke.c" function="cpapke_dec">
	  poly dec_vprime, dec_uhat, dec_tmp, dec_shat;

	  poly_frombytes(&dec_shat, sk);

	  decode_c(&dec_uhat, &dec_vprime, ct);
	  poly_mul_pointwise(&dec_tmp, &dec_shat, &dec_uhat);
	  poly_invntt(&dec_tmp);

	  poly_sub(&dec_tmp, &dec_tmp, &dec_vprime);

	  //poly_tomsg(out, &dec_tmp);
	  {
		// <inspired-by file="poly.c" function="poly_tomsg">
		unsigned int i;
		uint16_t t;

		for(i=0;i<32;i++)
		  out[i] = 0;

		for(i=0;i<256;i++)
		{
			if(i == target_coefficient)
				t  = flipabs_faulted((&dec_tmp)->coeffs[i+  0]);
			else
				t  = flipabs((&dec_tmp)->coeffs[i+  0]);

		  t += flipabs((&dec_tmp)->coeffs[i+256]);
	  #if (NEWHOPE_N == 1024)
		  t += flipabs((&dec_tmp)->coeffs[i+512]);
		  t += flipabs((&dec_tmp)->coeffs[i+768]);
		  t = ((t - NEWHOPE_Q));
	  #else
		  t = ((t - NEWHOPE_Q/2));
	  #endif

		  t >>= 15;
		  out[i>>3] |= t<<(i&7);
		}
		// </inspired-by>
	  }
	}

	int ineffective = (memcmp(&out, &msg, CRYPTO_BYTES) == 0);

//    fprintPoly(stdout, "compression_error = ", &compression_error);

	// write RE1_FILE
//	printf("\nr: ");
//	for (int i = 0; i < NEWHOPE_N; i++)
//	  printf("%d ", center(r.coeffs[i]));
//	printf("\ne1: ");
//	for (int i = 0; i < NEWHOPE_N; i++)
//	  printf("%d ", center(e_1.coeffs[i]));
//	printf("\n");
	print_linear_coeffs(ineqs_re1, &r, &e_1, target_coefficient);
	fprintf(ineqs_re1, " %s", ineffective ? ">=" : "<");
	fprintf(ineqs_re1, " %d\n", rhs);

	// print data
	/*fprintBstr(stdout, "msg = ", msg, CRYPTO_BYTES);
	fprintBstr(stdout, "ct = ", ct, CRYPTO_CIPHERTEXTBYTES);
	fprintBstr(stdout, "out = ", out, CRYPTO_BYTES);*/

	int16_t lhs = scalar_product(&r, &e_1, &e, &s, target_coefficient);


//	printf("%d %s %d\n", lhs, (ineffective ? ">=" : "<"), rhs);
	if (ineffective && (lhs - rhs < 0)) {
	  printf("FAIL assert(lhs >= rhs) failed\n");
	  return 3;
	} else if (!ineffective && (lhs - rhs >= 0)){
	  printf("FAIL assert(lhs < rhs) failed\n");
	  return 3;
	}
  }

  fclose(ineqs_es);
  fclose(ineqs_re1);
  return 0;
}


int main(int argc, char **argv){
	int numeqs = 10000;
	int coeff = 0;
	int filter = 100;

	int opt;

	// m: masked flag (default false)
	// s: seed (switches to deterministic mode)
	// n: number of equations
	// c: coefficient to attack
	// f: filter for ineqs, to only get those which have smaller than |f| on the rhs
	while((opt = getopt(argc, argv, "n:c:f:")) != -1){
		switch(opt){
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


	return generate_system(numeqs, coeff, filter);
}
