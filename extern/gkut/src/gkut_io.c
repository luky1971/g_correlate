/*
 * Copyright 2016 Ahnaf Siddiqui
 *
 * This program uses the GROMACS molecular simulation package API.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed at http://www.gromacs.org.
 */

#include "gkut_io.h"

#ifdef GRO_V5
#include "atoms.h"
#include "index.h"
#include "topology.h"
#include "trxio.h"
#endif
#include "smalloc.h"
#include "tpxio.h"
#include "typedefs.h"

void read_traj_t(const char *traj_fname, real **t, rvec ***x, matrix **box, int *nframes, int *natoms, output_env_t *oenv) {
	t_trxstatus *status = NULL;
	int est_frames = FRAMESTEP;
	*nframes = 0;

	snew(*t, est_frames);
	snew(*x, est_frames);
	snew(*box, est_frames);
	// *natoms = read_first_x(*oenv, &status, traj_fname, &((*t)[0]), &((*x)[0]), (*box)[0]);
	*natoms = read_first_x(*oenv, &status, traj_fname, *t, *x, (*box)[0]);

	do {
		++(*nframes);
		if(*nframes >= est_frames) {
			est_frames += FRAMESTEP;
			srenew(*t, est_frames);
			srenew(*x, est_frames);
			srenew(*box, est_frames);
		}
		snew((*x)[*nframes], *natoms);
	} while(read_next_x(*oenv, status, *t + *nframes,
#ifndef GRO_V5 
		*natoms,
#endif
		(*x)[*nframes], (*box)[*nframes]));

	sfree((*x)[*nframes]); // Nothing was read to the last allocated frame
	close_trx(status);
}

void read_traj(const char *traj_fname, rvec ***x, matrix **box, int *nframes, int *natoms, output_env_t *oenv) {
	real *t;
	read_traj_t(traj_fname, &t, x, box, nframes, natoms, oenv);
	sfree(t);
}

void print_traj(rvec **x, int nframes, int natoms, const char *fname) {
	int fr, i;
	FILE *f = fopen(fname, "w");

	for(fr = 0; fr < nframes; ++fr) {
		fprintf(f, "\nFrame %d\n", fr);
		for(i = 0; i < natoms; ++i) {
			fprintf(f, "%d: %f %f %f\n", i, x[fr][i][XX], x[fr][i][YY], x[fr][i][ZZ]);
		}
	}

	fclose(f);
}

void ndx_get_indx(const char *ndx_fname, int numgroups, atom_id ***indx, int **isize) {
	char **grp_names;

	snew(*isize, numgroups);
	snew(*indx, numgroups);
	snew(grp_names, numgroups);

	rd_index(ndx_fname, numgroups, *isize, *indx, grp_names);
	sfree(grp_names);
}

void filter_vecs(atom_id *indx, int isize, rvec *pre_x, rvec **new_x) {
	snew(*new_x, isize);
	for(int i = 0; i < isize; ++i) {
		copy_rvec(pre_x[indx[i]], (*new_x)[i]);
	}
}

void ndx_filter_traj(const char *ndx_fname, rvec **pre_x, rvec ***new_x, int nframes, int *natoms) {
	const int NUMGROUPS = 1;
	int *isize;
	atom_id **indx;

	ndx_get_indx(ndx_fname, NUMGROUPS, &indx, &isize);

	*natoms = isize[0];
	sfree(isize);

	snew(*new_x, nframes);
	for(int i = 0; i < nframes; ++i) {
		filter_vecs(indx[0], *natoms, pre_x[i], &((*new_x)[i]));
	}

	sfree(indx[0]);
	sfree(indx);
}


void read_top_gro(const char *gro_fname, t_topology *top) {
	char title[256];
	rvec *x = NULL;
	matrix box;
	int ePBC;

	init_top(top);

	read_tps_conf(gro_fname, title, top, &ePBC, &x, NULL, box, FALSE);

	sfree(x);
}

void free_topology(t_topology *top) {
	// Cannot use done_top(), causes error- pointer being freed was not allocated. See implementation in typedefs.c
	done_atom(&(top->atoms));
	done_symtab(&(top->symtab));
	done_block(&(top->cgs));
	done_block(&(top->mols));
	done_blocka(&(top->excls));
}

void read_top_tpr(const char *tpr_fname, gmx_mtop_t *mtop) {
	t_inputrec ir;
	matrix box;
	int natoms;

	read_tpx(tpr_fname, &ir, box, &natoms, NULL, NULL, NULL, mtop);
	done_inputrec(&ir);
}

void free_mtop(gmx_mtop_t *mtop) {
	done_mtop(mtop, TRUE);
}
