/*
 * rounding.cpp
 *
 *  Created on: Mar 16, 2020
 *      Author: ppessl
 */

#include <iostream>
#include <iomanip>
#include <cstdint>
#include <array>

extern "C"{
#include "kyber/api.h"
}


using namespace std;


void round_dist_v(int k=3){
	const auto v = 300;
	array<uint64_t, 2*v+1> hist = {0};

	for(auto i = 0; i < KYBER_Q; i++){
		int32_t error = 0;
		switch(k){
		case 2:{
			//relevant code from poly.c, line 28
			//t[j] = ((((uint32_t)a->coeffs[i+j] << 3) + KYBER_Q/2) / KYBER_Q) & 7
			auto compressed = ((((uint32_t)i << 3) + KYBER_Q/2) / KYBER_Q) & 7;

			//relevant code from poly.c, line 80
			//r->coeffs[i+0] =  (((a[0] & 7) * KYBER_Q) + 4) >> 3;
			auto decompressed = (((compressed & 7) * KYBER_Q) + 4) >> 3;

			error = (int32_t)decompressed - i;
			break;
		}
		case 3:{
			//relevant code from poly.c, line 39
			//t[j] = ((((uint32_t)a->coeffs[i+j] << 4) + KYBER_Q/2) / KYBER_Q) & 15;
			auto compressed = ((((uint32_t)i << 4) + KYBER_Q/2) / KYBER_Q) & 15;

			//relevant code from poly.c, line 93
			//r->coeffs[i+0] = (((a[0] & 15) * KYBER_Q) + 8) >> 4;
			auto decompressed = (((compressed & 15) * KYBER_Q) + 8) >> 4;

			error = (int32_t)decompressed - i;
			break;
		}
		case 4:{
			//relevant code from poly.c, line 51
		    //t[j] = ((((uint32_t)a->coeffs[i+j] << 5) + KYBER_Q/2) / KYBER_Q) & 31;
			auto compressed = ((((uint32_t)i << 5) + KYBER_Q/2) / KYBER_Q) & 31;

			//relevant code from poly.c, line 106
		    //r->coeffs[i+0] =  (((a[0] & 31) * KYBER_Q) + 16) >> 5;
			auto decompressed = (((compressed & 31) * KYBER_Q) + 16) >> 5;

			error = (int32_t)decompressed - i;
			break;
		}
		}


		if(error < -KYBER_Q/2)
			error += KYBER_Q;

		hist[error+v]++;

//		cout << setw(4) << i << " " << setw(4) << decompressed << endl;
	}

	for(auto i = 0; i < 2*v+1; i++){
		cout << setw(4) << i - v << " " << hist[i] << endl;
	}
}

void round_dist_u(int k=3){
	const auto v = 2;
	array<uint64_t, 2*v+1> hist = {0};

	for(auto i = 0; i < KYBER_Q; i++){
		int32_t error = 0;
		switch(k){
		case 2: case 3: {
			//relevant code from polyvec.c, line 49
			//t[k] = ((((uint32_t)a->vec[i].coeffs[4*j+k] << 10) + KYBER_Q/2) / KYBER_Q) & 0x3ff;
			auto compressed = ((((uint32_t)i << 10) + KYBER_Q/2) / KYBER_Q) & 0x3ff;

			//relevant code from polyvec.c, line 97
			//r->vec[i].coeffs[4*j+0] =  (((a[5*j+ 0]       | (((uint32_t)a[5*j+ 1] & 0x03) << 8)) * KYBER_Q) + 512) >> 10;
			auto decompressed = ((compressed * KYBER_Q) + 512) >> 10;

			error = (int32_t)decompressed - i;
			break;
		}
		case 4: {
			//relevant code from polyvec.c, line 26
			//t[k] = ((((uint32_t)a->vec[i].coeffs[8*j+k] << 11) + KYBER_Q/2) / KYBER_Q) & 0x7ff;
			auto compressed = ((((uint32_t)i << 11) + KYBER_Q/2) / KYBER_Q) & 0x7ff;

			//relevant code from polyvec.c, line 81
			//r->vec[i].coeffs[8*j+0] =  (((a[11*j+ 0]       | (((uint32_t)a[11*j+ 1] & 0x07) << 8)) * KYBER_Q) + 1024) >> 11;
			auto decompressed = ((compressed * KYBER_Q) + 1024) >> 11;

			error = (int32_t)decompressed - i;
			break;
		}
		}

		if(error < -KYBER_Q/2)
			error += KYBER_Q;

		hist[error+v]++;

//		cout << setw(4) << i << " " << setw(4) << decompressed << endl;
	}

	for(auto i = 0; i < 2*v+1; i++){
		cout << setw(4) << i - v << " " << hist[i] << endl;
	}
}
