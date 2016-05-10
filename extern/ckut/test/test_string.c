/*
 * Copyright 2016 Ahnaf Siddiqui
 */

#include <assert.h>
#include <stdio.h>
#include "ckut_string.h"

int main(int argc, char **argv) {
	assert(ck_strmatch("a", "a") == 0);
	assert(ck_strmatch("p", "a") != 0);
	assert(ck_strmatch("a", "p") != 0);
	assert(ck_strmatch("ahnaf", "ahnaf") == 0);
	assert(ck_strmatch("ab", "abc") != 0);
	assert(ck_strmatch("abc", "swag") != 0);
	assert(ck_strmatch("q*", "q") == 0);
	assert(ck_strmatch("Q*", "Quack") == 0);
	assert(ck_strmatch("Q*", "pewpew") != 0);
	assert(ck_strmatch("q*", "pq") != 0);
	assert(ck_strmatch("*q", "pq") == 0);
	assert(ck_strmatch("*o", "ayylmao") == 0);
	assert(ck_strmatch("N*2", "ND2") == 0);
	assert(ck_strmatch("N*2", "ND3") != 0);
	assert(ck_strmatch("N*2", "NDE2") == 0);

	printf("String tests passed.\n");
	return 0;
}
