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
#include "gkut_io.h"
#include "gkut_log.h"

#include "mtop_util.h"
#include "smalloc.h"

void get_corr_vecs(const char *top_fname) {
    t_topology top;

    switch(fn2ftp(top_fname)) {
        case efGRO:
                gk_read_top_gro(top_fname, &top);
            break;
        case efTPR:
            {
                gmx_mtop_t mtop;
                gk_read_top_tpr(top_fname, &mtop);
                top = gmx_mtop_t_to_t_topology(&mtop);
            }
            break;
        default:
            gk_log_fatal(FARGS, "%s is not a supported filetype for topology information!\n", top_fname);
    }

    if(top.idef.ntypes <= 0) {
        gk_free_topology(&top);
        gk_log_fatal(FARGS, "No interaction information found in %s!\n", top_fname);
    }

    // DEBUG
    for(int i = 0; i < top.atoms.nr; ++i) {
        printf("%d: %s-%d\n", i, top.atoms.atom[i].elem, top.atoms.atom[i].atomnumber);
    }

    // Free topology data
    gk_free_topology(&top);
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
        //  printf("%f\n", t[f]);
        // }

        // TODO: whole trajectory autocorrelation

        // Cleanup
        sfree(t);
        gk_free_traj(x, nframes, natoms);
        sfree(box);
    }
}
