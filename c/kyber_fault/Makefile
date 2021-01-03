SHELL=/bin/bash

KYBER_K = 2
KX = $$(( $(KYBER_K) * 256 ))

CC = gcc
CFLAGS = -Wall -Wextra -march=native -mtune=native -O3 -fomit-frame-pointer -DKYBER_K=$(KYBER_K)

CXX = g++
CXXFLAGS = $(CFLAGS) -std=c++0x

KSRCFOLDER = kyber
KSOURCES = aes256ctr.c cbd.c fips202.c indcpa.c kem.c kex.c ntt.c poly.c polyvec.c randombytes.c reduce.c sha256.c sha512.c symmetric-fips202.c verify.c
KHEADERS = aes256ctr.h api.h cbd.h fips202.h indcpa.h kex.h ntt.h params.h poly.h polyvec.h randombytes.h reduce.h sha2.h symmetric.h verify.h
KSOURCES := $(addprefix $(KSRCFOLDER)/, $(KSOURCES))
KHEADERS := $(addprefix $(KSRCFOLDER)/, $(KHEADERS))

SOURCES = basictests.cpp main.cpp maskeddec.cpp rounding.cpp
HEADERS = basictests.h maskeddec.h

KSOURCES += fixed_keys.c
KHEADERS += fixed_keys.h

.PHONY: .kyber kyber_fault cleano clean all
.DEFAULT_GOAL := all 

.kyber: $(KSOURCES) $(KHEADERS)
	$(CC) $(CFLAGS) $(KSOURCES) -std=c11 -c
	
kyber_fault: $(SOURCES) $(HEADERS) .kyber
	$(CXX) $(CXXFLAGS) *.o $(SOURCES) -o kyber$(KX)_fault
	
cleano:
	rm -f *.o
	
#hack, call make in make
all: 
	make kyber_fault KYBER_K=2
	make cleano
	make kyber_fault KYBER_K=3
	make cleano
	make kyber_fault KYBER_K=4
	make cleano

clean:
	rm -f kyber*_fault
	rm -f *.o
