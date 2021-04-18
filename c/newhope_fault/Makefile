CC = /usr/bin/gcc
LDFLAGS = 
CCFLAGS = -O3

#MAIN = PQCgenKAT_kem.c
MAIN_GEN = generate_lineqsys.c
SOURCES = kem.c cpapke.c fips202.c ntt.c poly.c precomp.c reduce.c rng.c verify.c
HEADERS = api.h cpapke.h fips202.h ntt.h params.h poly.h reduce.h rng.h verify.h

all: newhope512_fault

newhope512_fault: $(SOURCES) $(MAIN_GEN) $(HEADERS)
	$(CC) $(SOURCES) $(MAIN_GEN) -DNEWHOPE_N=512 $(CCFLAGS) -o $@


# needs more analysis
#newhope1024_fault: $(SOURCES) $(MAIN_GEN) $(HEADERS)
#	$(CC) $(SOURCES) $(MAIN_GEN) -DNEWHOPE_N 1024 $(CCFLAGS) -o $@

.PHONY: clean

clean:
	-rm -f newhope*_fault
