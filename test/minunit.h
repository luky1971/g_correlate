/*
 * These cute but useful macros were stolen from Jera Design (http://www.jera.com/techinfo/jtns/jtn002.html)
 * with modifications
 */

// Colors! (thank you Andrejs Cainikovs from http://stackoverflow.com/questions/3219393/stdlib-and-colored-output-in-c)
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define MU_PASS_COLOR	ANSI_COLOR_GREEN
#define MU_FAIL_COLOR	ANSI_COLOR_RED

#define PASS_COLOR


#define mu_assert(message, test) do { if (!(test)) return message; } while (0)

#define mu_run_test(test) do { char *message = test(); tests_run++; \
                                if (message) return message; } while (0)

extern int tests_run;