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
	// Limit the amount of trajectory frames loaded into memory.
	// When C_MEM_LIMIT is not set, the whole trajectory is loaded at once.
	C_MEM_LIMIT = 2,
};

struct res_corr_t {
	const char *res_name; // name of this residue
	real **auto_corr; // autocorrelation function values per atom, size [sum(natoms)][corr_dat.ncorr]
	real *s2; // S2 order parameter for each atom, size [sum(natoms)]
	int *natoms; // number of atoms in this residue for each atom type in corr_dat.atomtypes. Size corr_dat.ntypes.
	// The atoms in auto_corr[] and s2[] are grouped by atom type in the same order as they are specified in corr_dat.atomtypes.
	// natoms[] is also in the same order of atom types.
};

struct corr_dat_t {
	real *t; // time delays in autocorrelation function (domain) of size ncorr.
	struct res_corr_t *res_corr; // array of structs holding autocorrelation function values and S2 parameters for each residue. Size nres.
	// TODO: atomtypes and ntypes
	int nres; // number of residues
	int ncorr; // number of autocorrelation values per residue
};

void calc_ac(const char *fnames[], output_env_t *oenv, struct corr_dat_t *corr, unsigned long flags);
/* Calculates the autocorrelation funtions for the trajectory in fnames[efT_TRAJ] using the topology in fnames[efT_TOP].
 */

#endif // CORRELATE_H