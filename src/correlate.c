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
#include "ckut_string.h" // pattern matching atom names

#include "mtop_util.h" // dealing with topologies
#include "smalloc.h" // memory stuff
#include "vec.h" // vector ops for gc_get_unit_vecs()


#define GC_PAIR_ALLOC 10 // Initial amount by which to dynamically allocate array memory for atom pairs.
#define GC_FRAME_ALLOC 500 // Initial amount by which to dynamically allocate array memory of size # trajectory frames. 
#define GC_TIME_EPS 0.000001 // Epsilon for comparing floating point time values

// Test if two floating point values are nearly equivalent
#define GC_FLT_EQ(X, Y, EPS)    (((X) > ((Y) - EPS)) && ((X) < ((Y) + EPS)))
#define GC_TIME_EQ(X, Y)   GC_FLT_EQ(X, Y, GC_TIME_EPS)


void gc_init_corr_dat(struct gcorr_dat_t *corr) {
    corr->atomnames = NULL;
    corr->nnamepairs = 0;
    corr->dt = -1;
    corr->nt = -1;
    corr->found_atoms = NULL;
    corr->natompairs = NULL;
    corr->auto_corr = NULL;
    corr->s2 = NULL;
}

void gc_get_pairs(const t_atoms *atoms, const t_ilist *bonds,
               const char **atomnames, int nnamepairs,
               int natompairs[],
               int **pairs) {
    // Initialize number of pairs for each atomtype to zero.
    for(int i = 0; i < nnamepairs; ++i) {
        natompairs[i] = 0;
    }

    // Search for atom-atom pairs.
    int cur_alloc = GC_PAIR_ALLOC;
    int pair_ind = 0;
    snew(*pairs, cur_alloc);
    int pair_a, pair_b;
    // Loop through the given atom name pairs.
    for(int ai = 0; ai < nnamepairs; ++ai) {
        // Loop through all bond pairs in the system and search for instances of the current atom name pair.
        // This loop is run for each atom name pair to keep the found pairs in the same order as atomtypes.
        // This nested loop isn't as bad as it looks; the number of atom name pairs is usually only 1 or 2.
        for(int bi = 0; bi < bonds->nr; bi+=3) {
            pair_a = -1; // -1 means atom not found

            // Check if the current bond pair matches the current name pair in atomnames
            if(ck_strmatch(atomnames[2*ai], *(atoms->atomname[bonds->iatoms[bi+1]])) == 0 
                && ck_strmatch(atomnames[2*ai+1], *(atoms->atomname[bonds->iatoms[bi+2]])) == 0) {
                // DEBUG
                // printf("%s matches %s and %s matches %s\n", 
                //     atomnames[2*ai], *(atoms->atomname[bonds->iatoms[bi+1]]),
                //     atomnames[2*ai+1], *(atoms->atomname[bonds->iatoms[bi+2]]));

                pair_a = bonds->iatoms[bi+1];
                pair_b = bonds->iatoms[bi+2];
            }
            else if(ck_strmatch(atomnames[2*ai+1], *(atoms->atomname[bonds->iatoms[bi+1]])) == 0 
                && ck_strmatch(atomnames[2*ai], *(atoms->atomname[bonds->iatoms[bi+2]])) == 0) {
                // DEBUG
                // printf("%s matches %s and %s matches %s\n", 
                //     atomnames[2*ai], *(atoms->atomname[bonds->iatoms[bi+2]]),
                //     atomnames[2*ai+1], *(atoms->atomname[bonds->iatoms[bi+1]]));

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
                ++(natompairs[ai]);
            }
        }
    }
    
    srenew(*pairs, pair_ind * 2); // Free any excess memory in pairs.
}


void gc_get_unit_vecs(const rvec x[], const int pairs[], int npairs, rvec unit_vecs[]) {
    rvec temp;
    for(int p = 0; p < npairs; ++p) {
        rvec_sub(x[pairs[2 * p + 1]], x[pairs[2 * p]], temp);
        unitv(temp, unit_vecs[p]);
    }
}


void gc_correlate(const char *fnames[], output_env_t *oenv, struct gcorr_dat_t *corr, unsigned long flags) {
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
    snew(corr->natompairs, corr->nnamepairs);

    gc_get_pairs(&top.atoms, &bonds, corr->atomnames, corr->nnamepairs, corr->natompairs, &(corr->found_atoms));

    int npairs_tot = 0;
    for(int i = 0; i < corr->nnamepairs; ++i) {
        npairs_tot += corr->natompairs[i];
    }

    // DEBUG
    // Print found atom-atom pairs.
    for(int i = 0; i < npairs_tot * 2; i+=2) {
        printf("Pair %d: %d and %d, %s and %s, atomic #s %d and %d\n",
            i/2, corr->found_atoms[i], corr->found_atoms[i+1], 
            *(top.atoms.atomname[corr->found_atoms[i]]), *(top.atoms.atomname[corr->found_atoms[i+1]]), 
            top.atoms.atom[corr->found_atoms[i]].atomnumber, top.atoms.atom[corr->found_atoms[i+1]].atomnumber);
    }

    // Done with topology
    gk_free_topology(&top);


    // Read trajectory and calculate junk

    if(flags & C_MEM_LIMIT) {
        // Load a frame at a time and do calculations
    }
    else {
        // Store needed data from whole trajectory then do calculations
        rvec *x;
        t_trxstatus *status = NULL;
        int natoms = 0;
        real t = 0.0;
        matrix box;

        int nframes = 0;

        rvec **unit_vecs;
        int cur_alloc = GC_FRAME_ALLOC;
        snew(unit_vecs, cur_alloc);

        natoms = read_first_x(*oenv, &status, fnames[efT_TRAJ], &t, &x, box);

        if(natoms > 0) {
            real cur_dt = 0.0;
            real last_t = t;

            if(corr->dt < 0) {
                // User didn't provide dt, so use default:
                // Use trajectory timestep as dt and just store data for all trajectory frames

                real last_dt = corr->dt;

                do {
                    cur_dt = t - last_t;
                    last_t = t;

                    // Check for consistency of trajectory timestep 
                    // (this is necessary for proper autocorrelation with default trajectory timestep, 
                    // otherwise specify a timestep in corr->dt that is at least as large as the largest timestep
                    // in the trajectory and is a multiple of all other timesteps in the trajectory)
                    if(!GC_TIME_EQ(cur_dt, corr->dt) && nframes > 1) {
                        gk_log_fatal(FARGS, "Inconsistent time step in frame %d of %s: dt is %f vs. previous dt of %f.\n", 
                            nframes, fnames[efT_TRAJ], cur_dt, corr->dt);
                    }
                    else if(nframes == 1) {
                        // Use the dt between the first two frames of the trajectory as the expected dt.
                        corr->dt = cur_dt;
                    }

                    // Expand memory for unit vectors if needed
                    if(nframes + 1 > cur_alloc) {
                        cur_alloc *= 2;
                        srenew(unit_vecs, cur_alloc);
                    }

                    // Get unit vectors for the atom pairs in this frame
                    snew(unit_vecs[nframes], npairs_tot);

                    gc_get_unit_vecs(x, corr->found_atoms, npairs_tot, unit_vecs[nframes]);

                    ++nframes;

                } while(read_next_x(*oenv, status, &t,
#ifndef GRO_V5 
                    natoms,
#endif
                    x, box));

                if(corr->nt < 0) {
                    // User didn't provide nt, so use default which is
                    // maximum possible number of correlation intervals
                    corr->nt = nframes - 1;
                }

            } // if dt < 0
            else { // if provided dt >= 0
                // Only store the data in the trajectory intervals matching the given dt
                /*
                int cur_vec_ind = 0;

                do {
                    cur_dt += t - last_t;
                    last_t = t;

                    if(GC_TIME_EQ(cur_dt, corr->dt))

                    ++nframes;
                }
                */
            } // if provided dt >= 0
        } // if natoms > 0
        else {
            gk_log_fatal(FARGS, "No atoms found in %s!\n", fnames[efT_TRAJ]);
        }

        srenew(unit_vecs, nframes); // free excess memory

        // DEBUG
        gk_print_log("dt is %f, nt is %f.\n", corr->dt, corr->nt);
        // print unit vecs
        FILE *vecf = fopen("vecs.txt", "w");
        for(int fr = 0; fr < nframes; ++fr) {
            fprintf(vecf, "\nFrame %d:\n", fr);
            for(int p = 0; p < npairs_tot; ++p) {
                fprintf(vecf, "Pair %d: %f, %f, %f\n", p, unit_vecs[fr][p][XX], unit_vecs[fr][p][YY], unit_vecs[fr][p][ZZ]);
            }
        }
        fclose(vecf);
        gk_print_log("Unit vectors saved to vecs.txt for debugging.\n");

        // TODO: calculate autocorrelation from unit vectors



        // Cleanup

        sfree(x);

        for(int i = 0; i < nframes; ++i)
            sfree(unit_vecs[i]);
        sfree(unit_vecs);
    } // if not memory limit
}

// WARNING: this does not free any memory allocated for corr->atomnames.
// Whoever allocated that memory is responsible for it! 
void gc_free_corr(struct gcorr_dat_t *corr) {
    if(corr->found_atoms)   sfree(corr->found_atoms);

    if(corr->auto_corr) {
        int total = 0;
        for(int i = 0; i < corr->nnamepairs; ++i) {
            total += corr->natompairs[i];
        }
        for(int i = 0; i < total; ++i) {
            sfree(corr->auto_corr[i]);
        }
        sfree(corr->auto_corr);
    }

    if(corr->natompairs)    sfree(corr->natompairs);
    if(corr->s2)            sfree(corr->s2);
}
