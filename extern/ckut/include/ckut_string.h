/*
 * Copyright 2016 Ahnaf Siddiqui
 */

#ifndef CKUT_STRING_H
#define CKUT_STRING_H

#include <fnmatch.h>

#define FNMATCH_FLAGS 0

/**
  Checks whether the given string matches the given pattern.
  Returns 0 if the string matches the pattern, otherwise returns a nonzero value.
 */
static inline int ck_strmatch(const char *pattern, const char *string) {
	return fnmatch(pattern, string, FNMATCH_FLAGS);
}

#endif // CKUT_STRING_H