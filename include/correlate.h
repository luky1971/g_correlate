/*
 * Copyright 2016 Ahnaf Siddiqui and Sameer Varma
 *
 * This program uses the GROMACS molecular simulation package API.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed at http://www.gromacs.org.
 */

#ifndef CORRELATE_H
#define CORRELATE_H

#ifdef GRO_V5
#include "pargs.h"
#else
#include "statutil.h"
#endif

// Indices of filenames
enum {efT_TRAJ, efT_NDX, efT_TOP, efT_OUTDAT, efT_NUMFILES};

// Flags
enum { 
    C_FFT = 1, // Use FFT method for calculating autocorrelations.
    C_MEM_LIMIT = 2, // Limit the amount of trajectory frames loaded into memory.
                     // When C_MEM_LIMIT is not set, the whole trajectory is loaded at once.
};


struct corr_dat_t {
    const char **atomnames; // INPUT: the target atom name pairs to be tracked in autocorrelation. Size 2 * npairs.
                            // Should be formatted as a linear array of pairs, with two adjacent atom name strings corresponding to a pair.
                            // Supports wildcards.
                            // ex. given {"N", "H", "ND2", "H*", "C*", "H*"}, the program will search for N-H pairs, ND2-H* pairs, and C*-H* pairs,
                            // where H* will match any atom with a name starting with 'H', and likewise for C*.
    int npairs; // INPUT: the number of pairs in atomnames.

    real *t; // INPUT: time delays in autocorrelation function (domain), size ncorr.

    // The atom pairs in auto_corr[] and s2[] are grouped by atom name in the same order as they are specified in atomnames.
    // ie. for each atomnames pair i, there are natompairs[i] autocorrelation values,
    // followed by natompairs[i+1] autocorrelation values for pair i + 1, and so on until i + x = npairs.
    real **auto_corr; // autocorrelation function values for each atom pair, size [sum(natompairs)][ncorr]
    real *s2; // S2 order parameter for each atom pair, size [sum(natompairs)]
    int ncorr; // number of different time delays for autocorrelation.

    int *natompairs; // number of atom pairs found for each atom name pair in atomnames. Size npairs.

    int *res; // The residue ID of each atom pair, in same order as auto_corr and s2.
    const char **res_names; // The names of the residues. Indexed by the IDs in res. Size nres.
    int nres; // number of residues.
};


void calc_ac(const char *fnames[], output_env_t *oenv, struct corr_dat_t *corr, unsigned long flags);
/* Calculates the autocorrelation functions for the trajectory in fnames[efT_TRAJ] using the topology in fnames[efT_TOP].
 */

void free_corr(struct corr_dat_t *corr);
/* Frees the dynamic memory in a corr_dat_t struct.
 */

#endif // CORRELATE_H