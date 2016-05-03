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
    t_topology *top;

    /*
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
    }*/
    if(fn2ftp(top_fname) == efTOP || fn2ftp(top_fname) == efITP) {
        top = read_top(top_fname, NULL);
    }
    else {
        gk_log_fatal(FARGS, "%s is not a supported filetype for topology information!\n", top_fname);
    }

    if(top->idef.ntypes <= 0) {
        gk_free_topology(top);
        gk_log_fatal(FARGS, "No interaction information found in %s!\n", top_fname);
    }

    // DEBUG
    // Print atom names and atom numbers
    /*
    for(int i = 0; i < top.atoms.nr; ++i) {
        printf("%d: %s: %d\n", i, 
            *(top.atoms.atomname[i]), 
            top.atoms.atom[i].atomnumber);
    }
    */

    // DEBUG
    // Print counts for interaction types
    /*
    for(int f = 0; f < F_NRE; ++f) {
        printf("Interactions array for interaction type %d has %d elements.\n", f, top.idef.il[f].nr);
    }*/
    
    
    t_ilist bonds = top->idef.il[F_BONDS];
    printf("\nInteractions array for F_BONDS has %d elements.\n", bonds.nr);
    for(int i = 0; i < bonds.nr; i+=3) {
        printf("Bond %d, type %d: %d %d\n", i/3, bonds.iatoms[i], bonds.iatoms[i+1], bonds.iatoms[i+2]);
    }

    /*
    bonds = top.idef.il[F_G96BONDS];
    printf("\nInteractions array for G96BONDS has %d elements.\n", bonds.nr);
    for(int i = 0; i < bonds.nr; i+=3) {
        printf("Bond %d, type %d: %d %d\n", i/3, bonds.iatoms[i], bonds.iatoms[i+1], bonds.iatoms[i+2]);
    }*/
    
    // Free topology data
    gk_free_topology(top);
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
