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

// Command line parsing
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


struct gcorr_dat_t {
    const char **atomnames; // INPUT: the target atom name pairs to be tracked in autocorrelation. Size 2 * nnamepairs.
                            // Should be formatted as a linear array of pairs, with two adjacent atom name strings corresponding to a pair.
                            // Supports wildcards.
                            // ex. given {"N", "H", "ND2", "H*", "C*", "H*"}, the program will search for N-H pairs, ND2-H* pairs, and C*-H* pairs,
                            // where H* will match any atom with a name starting with 'H', and likewise for C*.
    int nnamepairs; // INPUT: the number of pairs in atomnames.
    real dt, nt; // INPUT: the time delay step and number of time delays in the domain of the autocorrelation functions.
                 // Set to -1 to use the default behavior. The default for dt is to use the trajectory's time step,
                 // and the default for nt is to go up to the length of the trajectory.

    int *found_atoms; // Gromacs IDs of the atoms found in the trajectory corresponding to the atom name pairs in atomnames.
                      // Size 2 * sum(natompairs). Members of a pair are adjacent.
                      // Atom pairs are grouped in the same order as the names in atomnames.
    int *natompairs; // number of atom pairs found for each atom name pair in atomnames, in same order. Size nnamepairs.

    // Each array in auto_corr[] and value in s2[] corresponds to an atom pair in the same order as found_atoms
    // and are grouped by atom name in the same order as they are specified in atomnames.
    // ie. for each atomnames pair i, there are natompairs[i] autocorrelation values in auto_corr[i],
    // followed by natompairs[i+1] autocorrelation values in auto_corr[i+1] for pair i + 1, and so on until i + x = sum(natompairs).
    real **auto_corr; // autocorrelation function values for each atom pair, size [sum(natompairs)][nt]
    real *s2; // S2 order parameter for each atom pair, size [sum(natompairs)]
};


void gc_init_corr_dat(struct gcorr_dat_t *corr);
/* Initializes a gcorr_dat_t struct, such as setting pointers to NULL and setting default parameters.
 */

void gc_get_pairs(const t_atoms *atoms, const t_ilist *bonds, // Input: Topology where atom-atom pairs will be searched.
               const char **atomnames, int nnamepairs, // Input: Pairs of atom names to be searched for in topology.
               int natompairs[], // Output: Number of atom pairs found for each atom name pair in atomnames.
                                 // Given buffer should be allocated to size nnamepairs.
               int **pairs); // Output 1D array of gromacs atom IDs of atom pairs. Memory is allocated for pairs.
                             // The pairs array will be size 2 * sum(natompairs).
                             // The members of a pair are adjacent.
                             // The pairs are grouped by atom names, in the same order as given atomnames.
                             // ie. for each atomnames pair i, there are natompairs[i] pairs of IDs in this pairs array,
                             // or 2 * natompairs[i] elements.
                             // The IDs in this array can be used to index into a gromacs trajectory associated with this topology.

void gc_get_unit_vecs(const rvec x[], const int pairs[], int npairs, rvec unit_vecs[]);
/* Calculates the unit vector formed between the points in each pair in the given array of points.
 * Pairs are specified by indexes into the array of points x.
 * The resulting unit vectors are stored in unit_vecs in the same order as pairs[].
 * unit_vecs should be pre-allocated to size npairs.
 */

void gc_correlate(const char *fnames[], output_env_t *oenv, struct gcorr_dat_t *corr, unsigned long flags);
/* Calculates the autocorrelation functions for the trajectory in fnames[efT_TRAJ] using the topology in fnames[efT_TOP].
 */

void gc_free_corr(struct gcorr_dat_t *corr);
/* Frees the dynamic memory in a gcorr_dat_t struct.
 */

#endif // CORRELATE_H