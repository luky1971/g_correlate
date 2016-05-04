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

#define ALLOC 10 // Initial amount by which to dynamically allocate array memory

void get_corr_pairs(const char *top_fname, // Topology file where atom-H pairs will be searched.
                    int atomtypes[], int n_atomtypes, // Atom types to be searched for in topology to find atom-H pairs.
                                                      // Identified by atomic number.
                    int natoms[], // Number of atoms of each type in atomtypes, in same order as atomtypes.
                                  // Given array should be size n_atomtypes.
                    int **pairs, int *npairs // Output array of atom IDs of atom-H pairs. Memory is allocated for pairs.
                                             // The first ID in a pair is of an input atomtype, the second is an associated hydrogen.
                                             // Therefore, the pairs array will be size 2 * npairs.
                                             // The pairs, grouped by the type of the first atom, 
                                             // are in the same order as the order of atomtypes.
                    ) {
    // Get topology data
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

    if(top.idef.ntypes <= 0 || ) {
        gk_free_topology(&top);
        gk_log_fatal(FARGS, "No interaction information found in %s!\n", top_fname);
    }

    t_ilist bonds = top.idef.il[F_BONDS];
    if(bonds.nr <= 0)
        bonds = top.idef.il[F_CONSTR];
    if(bonds.nr <= 0)
        gk_log_fatal(FARGS, "No bond information found in %s!\n", top_fname);


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
    
    // DEBUG
    // Print interactions
    /*
    printf("\nInteractions array has %d elements.\n", bonds.nr);
    for(int i = 0; i < bonds.nr; i+=3) {
        printf("Bond %d, type %d: %d %d\n", i/3, bonds.iatoms[i], bonds.iatoms[i+1], bonds.iatoms[i+2]);
    }*/

    // Search for atomtype-H pairs.
    int cur_alloc = ALLOC;
    
    for(int i = 0; i < n_atomtypes; ++i) {
        natoms[i] = 0;
    }

    int m_npairs = 0;
    int *mpairs = *pairs;
    snew(mpairs, cur_alloc);
    for(int ai = 0; ai < n_atomtypes; ++ai) {
        for(int bi = 0; bi < bonds.nr; bi+=3) {
            // if(bonds.iatoms[i+1])
        }
    }
    
    *npairs = m_npairs;

    // Free topology data
    gk_free_topology(&top);
}

void calc_ac(const char *fnames[], output_env_t *oenv, struct corr_dat_t *corr, unsigned long flags) {
    int *pairs, npairs;
    snew(corr->natoms, corr->n_atomtypes);

    get_corr_pairs(fnames[efT_TOP], corr->atomtypes, corr->n_atomtypes, corr->natoms, &pairs, &npairs);

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

    sfree(pairs);
    sfree(atomnums);
}

void free_corr(struct corr_dat_t *corr) {
    // TODO
}
