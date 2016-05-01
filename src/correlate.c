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
#include "gkut_log.h"

void get_corr_vecs(const char *top_fname) {
	switch(fn2ftp(top_fname)) {
		case efGRO:
			{
				t_topology top;

				gk_read_top_gro(top_fname, &top);
				/*
				if(top.idef.ntypes > 0) {
					ndih = top.idef.il[F_PDIHS].nr;

					read_traj(fnames[eTRAJ1], &x1, &nframes, &natoms, oenv);
					read_traj(fnames[eTRAJ2], &x2, &nframes, &natoms, oenv);

					ilist2svm_probs(top.idef.il, x1, x2, nframes, &probs);

					snew(models, ndih);
					train_traj(probs, ndih, gamma, c, parallel, models);

					// Construct eta_dih struct
					eta_dih->ndih = ndih;

					snew(eta_dih->atoms, ndih * 4);
					int i;
					for(i = 0; i < ndih; ++i) {
						eta_dih->atoms[4*i] = top.idef.il[F_PDIHS].iatoms[i*5+1];
						eta_dih->atoms[4*i+1] = top.idef.il[F_PDIHS].iatoms[i*5+2];
						eta_dih->atoms[4*i+2] = top.idef.il[F_PDIHS].iatoms[i*5+3];
						eta_dih->atoms[4*i+3] = top.idef.il[F_PDIHS].iatoms[i*5+4];
					}
				}
				else {
					print_log("No interaction information found in %s. Skipping dihedral eta.\n", fnames[eTOP1]);
					success = FALSE;
				}*/

				gk_free_topology(&top);
			}
			break;
		case efTPR:
			{
				gmx_mtop_t mtop;

				gk_read_top_tpr(top_fname, &mtop);
				/*
				ndih = mtop.moltype->ilist[F_PDIHS].nr;

				read_traj(fnames[eTRAJ1], &x1, &nframes, &natoms, oenv);
				read_traj(fnames[eTRAJ2], &x2, &nframes, &natoms, oenv);

				ilist2svm_probs(mtop.moltype->ilist, x1, x2, nframes, &probs);

				snew(models, ndih);
				train_traj(probs, ndih, gamma, c, parallel, models);

				// Construct eta_dih struct
				eta_dih->ndih = ndih;

				snew(eta_dih->atoms, ndih * 4);
				int i;
				for(i = 0; i < ndih; ++i) {
					eta_dih->atoms[4*i] = mtop.moltype->ilist[F_PDIHS].iatoms[i*5+1];
					eta_dih->atoms[4*i+1] = mtop.moltype->ilist[F_PDIHS].iatoms[i*5+2];
					eta_dih->atoms[4*i+2] = mtop.moltype->ilist[F_PDIHS].iatoms[i*5+3];
					eta_dih->atoms[4*i+3] = mtop.moltype->ilist[F_PDIHS].iatoms[i*5+4];
				}*/

				gk_free_mtop(&mtop);
			}
			break;
		default:
			gk_log_fatal(FARGS, "%s is not a supported filetype for topology information!\n", top_fname);
	}
}

void calc_ac(const char *fnames[], output_env_t *oenv, struct corr_dat_t *corr, unsigned long flags) {
	get_corr_vecs(fnames[efT_TOP]);

	if(flags & C_MEM_LIMIT) {
		// Load a frame at a time and do calculations
	}
	else {
		// Get whole trajectory then do calculations
		real *t;
		rvec **x;
		matrix *box;
		int nframes, natoms;

		gk_read_traj_t(fnames[efT_TRAJ], &t, &x, &box, &nframes, &natoms, oenv);

		// DEBUG
		// for(int f = 0; f < nframes; ++f) {
		// 	printf("%f\n", t[f]);
		// }

		// TODO: whole trajectory autocorrelation

		// Cleanup
		sfree(t);
		gk_free_traj(x, nframes, natoms);
	    sfree(box);
	}
}
