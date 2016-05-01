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

#include "correlate.h"
#include "smalloc.h"
#include "gkut_io.h"

void calc_ac(const char *fnames[], output_env_t *oenv, struct corr_dat_t *corr, unsigned long flags) {
	if(flags & C_MEM_LIMIT) {
		// Load a frame at a time and do calculations
	}
	else {
		// Get whole trajectory then do calculations
		real *t;
		rvec **x;
		matrix *box;
		int nframes, natoms;

		read_traj_t(fnames[efT_TRAJ], &t, &x, &box, &nframes, &natoms, oenv);

		// DEBUG
		// for(int f = 0; f < nframes; ++f) {
		// 	printf("%f\n", t[f]);
		// }

		// TODO: whole trajectory autocorrelation

		// Cleanup
		sfree(t);
		for(int i = 0; i < nframes; ++i) {
	        sfree(x[i]);
	    }
	    sfree(x);
	    sfree(box);
	}
}
