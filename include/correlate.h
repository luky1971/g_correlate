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


// The atoms in auto_corr[] and s2[] are grouped by atom type in the same order as they are specified in atomtypes.
// natoms[] is also in the same order of atom types.
struct corr_dat_t {
    int *atomtypes; // INPUT: the target atom types, identified by atomic number, to be tracked in autocorrelation.
    int n_atomtypes; // INPUT: the number of elements in atomtypes.

    real *t; // INPUT: time delays in autocorrelation function (domain), size ncorr.
    real **auto_corr; // autocorrelation function values for each atom-H pair, size [sum(natoms)][ncorr]
    real *s2; // S2 order parameter for each atom-H pair, size [sum(natoms)]
    int ncorr; // number of different time delays for autocorrelation.

    int *natoms; // number of atoms found for each atom type in atomtypes. Size n_atomtypes.

    int *res; // The residue ID of each atom, in same order as auto_corr and s2.
    const char **res_names; // The names of the residues indexed by the IDs in res.
    int nres; // number of residues.
};


void calc_ac(const char *fnames[], output_env_t *oenv, struct corr_dat_t *corr, unsigned long flags);
/* Calculates the autocorrelation funtions for the trajectory in fnames[efT_TRAJ] using the topology in fnames[efT_TOP].
 */

void free_corr(struct corr_dat_t *corr);
/* Frees the dynamic memory in a corr_dat_t struct.
 */

#endif // CORRELATE_H