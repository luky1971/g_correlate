
#include <stdio.h>
#include "correlate.h"
#include "minunit.h"

// Test counter
int tests_run = 0;

// Floating point equality tests
#define VEC_EPS 0.001

#define GC_FLT_EQ(X, Y, EPS)    (((X) > ((Y) - EPS)) && ((X) < ((Y) + EPS)))
#define GC_VEC_EQ(A, B)			(GC_FLT_EQ(A[0], B[0], VEC_EPS) && GC_FLT_EQ(A[1], B[1], VEC_EPS) && GC_FLT_EQ(A[2], B[2], VEC_EPS))


static char *test_unit_vecs() {
	rvec x[5];

	x[0][XX] = 1;
	x[0][YY] = 2;
	x[0][ZZ] = -4;
	x[1][XX] = -9;
	x[1][YY] = -1;
	x[1][ZZ] = 7;
	x[2][XX] = 15;
	x[2][YY] = 12;
	x[2][ZZ] = 6;
	x[3][XX] = -19;
	x[3][YY] = -1;
	x[3][ZZ] = -9;
	x[4][XX] = 3;
	x[4][YY] = -16;
	x[4][ZZ] = 0;

	int pairs[] = { 0, 2, 4, 1, 3, 4 };
#define NPAIRS 3

	rvec unit_vecs[NPAIRS];

	for(int p = 0; p < NPAIRS; ++p)
		gc_get_unit_vec(x, pairs[2*p], pairs[2*p+1], unit_vecs[p]);

	rvec exp[NPAIRS];

	exp[0][XX] = 0.703526;
	exp[0][YY] = 0.502519;
	exp[0][ZZ] = 0.502519;
	exp[1][XX] = -0.586939;
	exp[1][YY] = 0.733674;
	exp[1][ZZ] = 0.342381;
	exp[2][XX] = 0.782725;
	exp[2][YY] = -0.533676;
	exp[2][ZZ] = 0.320206;

	for(int i = 0; i < NPAIRS; ++i) {
		printf("Unit vector %d is %f, %f, %f; expected %f, %f, %f.\n", 
			i, unit_vecs[i][XX], unit_vecs[i][YY], unit_vecs[i][ZZ], exp[i][XX], exp[i][YY], exp[i][ZZ]);
	}

	for(int i = 0; i < NPAIRS; ++i) {
		printf("Testing pair %d...\n", i);
		mu_assert("bad unit vector!", GC_VEC_EQ(unit_vecs[i], exp[i]));
	}

	return NULL;
}

static char *all_tests() {
	mu_run_test(test_unit_vecs);
	return NULL;
}

int main(int argc, char **argv) {
	char *result = all_tests();
	if(result) {
		printf(MU_FAIL_COLOR "TEST FAILED" ANSI_COLOR_RESET ", %s\n", result);
	}
	else {
		printf(MU_PASS_COLOR "ALL TESTS PASSED!" ANSI_COLOR_RESET "\n");
	}
	printf("%d tests run.\n", tests_run);

	return result != NULL;
}
