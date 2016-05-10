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

#define GS2_ALLOC 10 // Initial amount by which to dynamically allocate array memory.

void get_corr_pairs(t_atoms *atoms, t_ilist *bonds, // Input: Topology where atom-atom pairs will be searched.
                    const char **atomnames, int npairs, // Input: Pairs of atom names to be searched for in topology.
                    int natompairs[], // Output: Number of atom pairs found for each atom name pair in atomnames.
                                      // Given buffer should be allocated to size npairs.
                    int **pairs // Output 1D array of gromacs atom IDs of atom pairs. Memory is allocated for pairs.
                                // The pairs array will be size 2 * sum(natompairs).
                                // The members of a pair are adjacent.
                                // The pairs are grouped by atom names, in the same order as given atomnames.
                                // ie. for each atomnames pair i, there are natompairs[i] pairs of IDs in this pairs array,
                                // or 2 * natompairs[i] elements.
                                // The IDs in this array can be used to index into a gromacs trajectory associated with this topology.
                    ) {
    // Initialize number of pairs for each atomtype to zero.
    for(int i = 0; i < npairs; ++i) {
        natoms[i] = 0;
    }

    // Search for atom-atom pairs.
    int cur_alloc = GS2_ALLOC;
    int pair_ind = 0;
    snew(*pairs, cur_alloc);
    int pair_a, pair_b;
    // Loop through the given atom name pairs.
    for(int ai = 0; ai < npairs; ++ai) {
        // Loop through all bond pairs in the system and search for instances of the current atom name pair.
        // This loop is run for each atom name pair to keep the found pairs in the same order as atomtypes.
        // This nested loop isn't as bad as it looks; the number of atom name pairs is usually only 1 or 2.
        for(int bi = 0; bi < bonds->nr; bi+=3) {
            pair_a = -1; // -1 means atom not found

            // Check if the current bond pair matches the current name pair in atomnames
            if(*(atoms->atomname[bonds->iatoms[bi+1]]) == atomnames[2*ai] 
                && *(atoms->atomname[bonds->iatoms[bi+2]]) == atomnames[2*ai+1]) {
                pair_a = bonds->iatoms[bi+1];
                pair_b = bonds->iatoms[bi+2];
            }
            else if(*(atoms->atomname[bonds->iatoms[bi+1]]) == atomnames[2*ai+1] 
                && *(atoms->atomname[bonds->iatoms[bi+2]]) == atomnames[2*ai]) {
                pair_a = bonds->iatoms[bi+2];
                pair_b = bonds->iatoms[bi+1];
            }

            // If matched, add the atom IDs of the current bond pair to pairs.
            if(pair_a >= 0) {
                // Found an atom-atom pair.
                // Increase size of pairs array if necessary.
                if(pair_ind > cur_alloc - 2) {
                    cur_alloc *= 2;
                    srenew(*pairs, cur_alloc);
                }
                // Add the found atom-atom pair to the pairs array.
                (*pairs)[pair_ind++] = pair_a;
                (*pairs)[pair_ind++] = pair_b;
                // Increment number of pairs for current atomtype.
                ++(natoms[ai]);
            }
        }
    }
    
    srenew(*pairs, pair_ind * 2); // Free any excess memory in pairs.
}


void calc_ac(const char *fnames[], output_env_t *oenv, struct corr_dat_t *corr, unsigned long flags) {
    // Get topology data
    t_topology top;

    switch(fn2ftp(fnames[efT_TOP])) {
        case efGRO:
                gk_read_top_gro(fnames[efT_TOP], &top);
            break;
        case efTPR:
            {
                gmx_mtop_t mtop;
                gk_read_top_tpr(fnames[efT_TOP], &mtop);
                top = gmx_mtop_t_to_t_topology(&mtop);
            }
            break;
        default:
            gk_log_fatal(FARGS, "%s is not a supported filetype for topology information!\n", fnames[efT_TOP]);
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
    
    // DEBUG
    // Print interactions
    /*
    printf("\nInteractions array has %d elements.\n", bonds.nr);
    for(int i = 0; i < bonds.nr; i+=3) {
        printf("Bond %d, type %d: %d %d\n", i/3, bonds.iatoms[i], bonds.iatoms[i+1], bonds.iatoms[i+2]);
    }*/

    if(top.idef.ntypes <= 0) {
        gk_free_topology(&top);
        gk_log_fatal(FARGS, "No interaction information found in %s!\n", fnames[efT_TOP]);
    }

    t_ilist bonds = top.idef.il[F_BONDS];
    if(bonds.nr <= 0)
        bonds = top.idef.il[F_CONSTR];
    if(bonds.nr <= 0) {
        gk_free_topology(&top);
        gk_log_fatal(FARGS, "No bond information found in %s!\n", fnames[efT_TOP]);
    }


    // Get atom-atom pairs
    snew(corr->natoms, corr->npairs);

    get_corr_pairs(&top.atoms, &bonds, corr->atomtypes, corr->npairs, corr->natoms, &(corr->found_atoms));

    int total = 0;
    for(int i = 0; i < corr->npairs; ++i) {
        total += corr->natoms[i];
    }

    // DEBUG
    // Print found atom-atom pairs.
    for(int i = 0; i < total * 2; i+=2) {
        printf("Pair %d: %d and %d, %s and %s, atomic #s %d and %d\n",
            i/2, corr->found_atoms[i], corr->found_atoms[i+1], 
            *(top.atoms.atomname[corr->found_atoms[i]]), *(top.atoms.atomname[corr->found_atoms[i+1]]), 
            top.atoms.atom[corr->found_atoms[i]].atomnumber, top.atoms.atom[corr->found_atoms[i+1]].atomnumber);
    }

    // Read trajectory and calculate junk
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

    // Cleanup
    gk_free_topology(&top);
}

void free_corr(struct corr_dat_t *corr) {
    // TODO
}
