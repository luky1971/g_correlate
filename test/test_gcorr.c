
#include <stdio.h>
#include "minunit.h"

int tests_run = 0;

static char *test_unit_vecs() {
	mu_assert("true is false!", 1);
	return NULL;
}

static char *all_tests() {
	mu_run_test(test_unit_vecs);
	return NULL;
}

int main(int argc, char **argv) {
	char *result = all_tests();
	if(result) {
		printf("ERROR, %s\n", result);
	}
	else {
		printf("ALL TESTS PASSED!\n");
	}
	printf("%d tests run.\n", tests_run);

	return result != NULL;
}
