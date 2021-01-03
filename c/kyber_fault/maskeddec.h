/*
 * maskeddec.h
 *
 *  Created on: Mar 18, 2020
 *      Author: ppessl
 */

#ifndef MASKEDDEC_H_
#define MASKEDDEC_H_

#define QBITS (12)
#define QBITS2 (QBITS+1)
#define QBYTES (2)
#define MSB(x) ((x >> (QBITS2-1)) & 0x1)

#include <cstdint>

extern "C"{
#include "kyber/api.h"
#include "kyber/indcpa.h" //poly and polyvec
}

void mask_poly(poly& s1, poly& s2, const poly& s);
void mask_polyvec(polyvec& s1, polyvec& s2, const polyvec& s);
void refresh_poly(poly& s1, poly& s2);
void refresh_polyvec(polyvec& s1, polyvec& s2);
void masked_decode(uint16_t& m1, uint16_t& m2, int16_t c1, int16_t c2);

void poly_negate(poly* r);
void polyvec_negate(polyvec* r);
void test_masked();

void indcpa_dec_fault_masked(unsigned char *m, const unsigned char *c, polyvec* s1, polyvec* s2, int faultcoeff);


#endif /* MASKEDDEC_H_ */
