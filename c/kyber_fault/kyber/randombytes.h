#ifndef RANDOMBYTES_H
#define RANDOMBYTES_H

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include <unistd.h>

void randombytes(unsigned char *x, size_t xlen);

#endif
