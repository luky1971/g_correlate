/*
 * Copyright 2016 Ahnaf Siddiqui and Sameer Varma
 *
 * 
 * Autocorrelation and S2 order parameter formulas are borrowed from
 *
 * Chatfield, D. C.; Szabo, A.; Brooks, B. R., 
 * Molecular Dynamics of Staphylococcal Nuclease:  Comparison of Simulation with 15N and 13C NMR Relaxation Data. 
 * Journal of the American Chemical Society 1998, 120 (21), 5301-5311.
 *
 * and
 *
 * Gu, Y.; Li, D.-W.; Brüschweiler, R., 
 * NMR Order Parameter Determination from Long Molecular Dynamics Trajectories for Objective Comparison with Experiment. 
 * Journal of Chemical Theory and Computation 2014, 10 (6), 2599-2607.
 *
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

#include <ctype.h> // for isspace
#ifdef GRO_V5
#include "pargs.h"
#else
#include "statutil.h"
#endif // Command line parsing
#include "vec.h" // vector ops for gc_get_unit_vec()

// Indices of filenames
enum {efT_TRAJ, efT_NDX, efT_TOP, efT_OUTDAT, efT_S2DAT, efT_NUMFILES};

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
    real dt; // INPUT: the time delay step (dt) and number of time delays (nt) in the domain of the autocorrelation functions.
    int nt;  // Set either or both dt and nt to -1 to use the default behavior. The default for dt is to use the trajectory's time step,
             // and the default for nt is to go up to the length of the trajectory.
             // See auto_corr below for more info on dt and nt.

    int *atompairs; // Gromacs IDs of the atoms found in the trajectory corresponding to the atom name pairs in atomnames.
                    // Size 2 * sum(natompairs). Members of a pair are adjacent.
                    // Atom pairs are grouped in the same order as the names in atomnames.
    int *natompairs; // number of atom pairs found for each atom name pair in atomnames, in same order. Size nnamepairs.

    // Each unit vector in unit_vecs 
    // and each array in auto_corr[] and value in s2[] corresponds to an atom pair in the same order as atompairs
    // and are grouped by atom name in the same order as they are specified in atomnames.
    // ie. for each atomnames pair i, there are natompairs[i] unit vectors in unit_vecs, 
    // and there are natompairs[i] autocorrelation functions in auto_corr,
    // followed by natompairs[i+1] autocorrelation fucntions for pair i + 1, and so on until i + x = sum(natompairs).
    rvec **unit_vecs; // the bond orientation unit vector for each pair in atompairs. Size [sum(natompairs)][nframes]
    int nframes; // The number of time points at which unit vectors were recorded.
                 // This is less than or equal to the number of frames in the trajectory (see gc_traj2uvecs function)
    real **auto_corr; // autocorrelation function values for each atom pair, size [sum(natompairs)][nt]
                      // for each atom pair, the first autocorrelation value corresponds to t = dt, i.e. t = 0 is skipped,
                      // and each subsequent t is the previous t plus dt,
                      // so that the last autocorrelation value corresponds to t = dt * nt.
                      // The skipped autocorrelation value at t = 0 is implied to be 1,
                      // this is always the case because unit vectors are used.
    real *s2; // S2 order parameter for each atom pair, size [sum(natompairs)]
};


void gc_init_corr_dat(struct gcorr_dat_t *corr);
/* Initializes a gcorr_dat_t struct, such as setting pointers to NULL and setting default parameters.
 */

void gc_correlate(const char *fnames[], output_env_t *oenv, struct gcorr_dat_t *corr, unsigned long flags);
/* Calculates the autocorrelation functions for the trajectory in fnames[efT_TRAJ] using the topology in fnames[efT_TOP].
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

int gc_traj2uvecs(const char *traj_fname, 
                  output_env_t *oenv, 
                  real *dt, 
                  int atompairs[], int npairs, 
                  rvec ***unit_vecs);
/* Calculates the unit vectors between the pairs of atoms given in atompairs in the trajectory traj_fname,
 * storing the resulting vectors in unit_vecs and returning the number of frames of vectors stored.
 * Memory is allocated for unit_vecs, which will be size nframes x npairs, where nframes is the function's return value.
 * The number of frames returned is not necessarily equal to the number of frames in the trajectory, 
 * as the frames between every dt timestep are skipped.
 * dt should NOT be NULL. The pointed value should be -1 to use the trajectory's timestep, which will be saved to dt.
 * Otherwise, the given dt value should be a multiple of the trajectory's timestep.
 */

void gc_calc_ac(const rvec vecs[], int nvecs, int nt, real auto_corr[]);
/* Calculates the autocorrelation function for a trajectory of vectors
 * and stores the results in auto_corr. auto_corr should be pre-allocated to size nt.
 */

real gc_calc_s2(const rvec unit_vecs[], int nvecs);
/* Calculates the generalized Lipari-Szabo order parameter S^2 for the given trajectory of unit bond vectors.
 * Returns S^2.
 */

void gc_save_corr(struct gcorr_dat_t *corr, const char *corr_fname, const char *s2_fname);
/* Outputs the given autocorrelation and s2 data to files with the given names.
 * Either corr_fname or s2_fname can be NULL to not output the corresponding data.
 */

void gc_free_corr(struct gcorr_dat_t *corr);
/* Frees the dynamic memory in a gcorr_dat_t struct.
 * WARNING: this does not free any memory allocated for corr->atomnames.
 * Whoever allocated that memory is responsible for it! 
 */


/* Calculates the unit vector formed between the points given by IDs a and b in the given array of points x.
 * a and b are indexes into the array of points x.
 * The resulting unit vector is stored in unit_vec.
 */
static inline void gc_get_unit_vec(const rvec x[], const int a, const int b, rvec unit_vec) {
    rvec temp;
    rvec_sub(x[b], x[a], temp);
    unitv(temp, unit_vec);
}


/* Parses the inwords string inline and stores pointers to parsed substrings in outwords, which must be pre-allocated.
 */
// NOTE: this modifies the inwords array and simply sets pointers to different portions of it in outwords,
// to avoid allocating and copying memory.
static void inparse(char *inwords, const char **outwords, int *nwords) {
    char *c = inwords;
    const char **w = outwords;
    int newword = 1;
    while(c != NULL && *c != '\0') {
        if(isspace(*c) || *c == '-' || *c == ',') {
            *c = '\0';
            newword = 1;
        }
        else if(newword) {
            *(w++) = c;
            newword = 0;
        }
        ++c;
    }
    *nwords = w - outwords;
}

#endif // CORRELATE_H