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

#include "correlate.h"
#include "gkut_io.h"
#include "gkut_log.h"
#include "ckut_string.h" // pattern matching atom names

#include "mtop_util.h" // dealing with topologies
#include "smalloc.h" // memory stuff


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
    corr->atompairs = NULL;
    corr->natompairs = NULL;
    corr->auto_corr = NULL;
    corr->s2 = NULL;
}


void gc_get_pairs(const t_atoms *atoms, const t_ilist *bonds,
               const char **atomnames, int nnamepairs,
               int natompairs[],
               int **pairs) {
    // Initialize number of pairs for each pair of atom names to zero.
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
        // This loop is run for each atom name pair to keep the found pairs in the same order as atomnames.
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
                // Increment number of pairs for current atomnames.
                ++(natompairs[ai]);
            }
        }
    }
    
    srenew(*pairs, pair_ind * 2); // Free any excess memory in pairs.
}


int gc_traj2uvecs(const char *traj_fname, 
				  output_env_t *oenv, 
				  real *dt, 
				  int atompairs[], int npairs, 
				  rvec ***unit_vecs) {
	rvec *x;
    t_trxstatus *status = NULL;
    int natoms = 0;
    real t = 0.0;
    matrix box;

    int nframes = 0;

    int cur_alloc = GC_FRAME_ALLOC;

    // allocate memory for unit vectors
    snew(*unit_vecs, cur_alloc);

    natoms = read_first_x(*oenv, &status, traj_fname, &t, &x, box);

    if(natoms > 0) {
        real cur_dt = 0.0;
        real last_t = t;

        if(*dt < 0) {
            // User didn't provide dt, so use default:
            // Use trajectory timestep as dt and just store data for all trajectory frames

            real last_dt = *dt;

            do {
                cur_dt = t - last_t;
                last_t = t;

                // Check for consistency of trajectory timestep 
                // (this is necessary for proper autocorrelation with default trajectory timestep, 
                // otherwise specify a timestep in dt that is at least as large as the largest timestep
                // in the trajectory and is a multiple of all other timesteps in the trajectory)
                if(!GC_TIME_EQ(cur_dt, *dt) && nframes > 1) {
                    gk_log_fatal(FARGS, "Inconsistent time step in frame %d of %s: dt is %f vs. previous dt of %f.\n", 
                        nframes, traj_fname, cur_dt, *dt);
                }
                else if(nframes == 1) {
                    // Use the dt between the first two frames of the trajectory as the expected dt.
                    *dt = cur_dt;
                }

                // Expand memory if needed
                if(nframes + 1 > cur_alloc) {
                	cur_alloc *= 2;
                	srenew(*unit_vecs, cur_alloc);
                }

                // Get unit vectors for the atom pairs in this frame
                snew((*unit_vecs)[nframes], npairs);
                for(int p = 0; p < npairs; ++p) {
                    gc_get_unit_vec(x, atompairs[2 * p], atompairs[2 * p + 1], (*unit_vecs)[nframes][p]);
                }

                ++nframes;

            } while(read_next_x(*oenv, status, &t,
#ifndef GRO_V5 
                natoms,
#endif
                x, box));
        } // if dt < 0
        else { // if provided dt >= 0
            // Only store the data in the trajectory intervals matching the given dt
            /*
            int cur_vec_ind = 0;

            do {
                cur_dt += t - last_t;
                last_t = t;

                if(GC_TIME_EQ(cur_dt, *dt))

                ++nframes;
            }
            */
        } // if provided dt >= 0
    } // if natoms > 0
    else {
        gk_log_fatal(FARGS, "No atoms found in %s!\n", traj_fname);
    }

    // Done with trajectory
    close_trx(status);
    sfree(x);

    srenew(*unit_vecs, nframes); // free excess memory

    return nframes;
}


void gc_calc_ac(const rvec *vecs, int nvecs, real nt, real *auto_corr) {
	// tdelay is the number of indexes of separation between consecutive vectors in the current time delay
	for(int tdelay = 1; tdelay <= nt; ++tdelay) { // iterate through time delays of autocorrelation domain
		real sum = 0, dot;
		for(int i = 0; i < nvecs - tdelay; i += tdelay) { // iterate through the vectors for each time delay
			dot = iprod(vecs[i], vecs[i + tdelay]);
			sum += 3.0 * dot * dot - 1.0;
		}
		auto_corr[tdelay - 1] = 0.5 * (sum / ((nvecs - 1) / tdelay));
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
    gc_get_pairs(&top.atoms, &bonds, corr->atomnames, corr->nnamepairs, corr->natompairs, &(corr->atompairs));


    // Calculate total number of pairs
    int npairs_tot = 0;
    for(int i = 0; i < corr->nnamepairs; ++i) {
        npairs_tot += corr->natompairs[i];
    }

    // DEBUG
    // Print found atom-atom pairs.
    for(int i = 0; i < npairs_tot * 2; i+=2) {
        printf("Pair %d: %d and %d, %s and %s, atomic #s %d and %d\n",
            i/2, corr->atompairs[i], corr->atompairs[i+1], 
            *(top.atoms.atomname[corr->atompairs[i]]), *(top.atoms.atomname[corr->atompairs[i+1]]), 
            top.atoms.atom[corr->atompairs[i]].atomnumber, top.atoms.atom[corr->atompairs[i+1]].atomnumber);
    }


    // Done with topology
    gk_free_topology(&top);


    // Read trajectory and get unit vectors for atom pairs
    rvec **unit_vecs;
    int nframes = gc_traj2uvecs(fnames[efT_TRAJ], oenv, &(corr->dt), corr->atompairs, npairs_tot, &unit_vecs);
    
    // Set default nt if needed
    if(corr->nt < 0) {
        // User didn't provide nt, so use default which is
        // maximum possible number of correlation intervals
        corr->nt = nframes - 1;
    }

    // DEBUG
    gk_print_log("dt is %f, nt is %d.\n", corr->dt, corr->nt);
    // print unit vecs
    FILE *vecf = fopen("vecs_fbyp.txt", "w");
    for(int fr = 0; fr < nframes; ++fr) {
        fprintf(vecf, "\nFrame %d:\n", fr);
        for(int p = 0; p < npairs_tot; ++p) {
            fprintf(vecf, "Pair %d: %f, %f, %f\n", p, unit_vecs[fr][p][XX], unit_vecs[fr][p][YY], unit_vecs[fr][p][ZZ]);
        }
    }
    fclose(vecf);
    gk_print_log("Unit vectors saved to vecs_fbyp.txt for debugging.\n");


    // Reformat 2d array unit_vecs from [frame#][vec] to [vec][frame#]
    // for better cache locality during autocorrelation calculations, 
    // since autocorrelation is calculated vector by vector.
    // Why wasn't unit_vecs formatted this way in the first place?
    // Because the reallocation pattern needed for growing the number of frames
    // in unit_vecs when reading the trajectory was causing memory issues.
    // This may be explored more in the future.

    // Create new 2d array for reformatted unit vectors
    rvec **old_unit_vecs = unit_vecs;
    snew(unit_vecs, npairs_tot);

    // Transpose unit_vecs matrix so that new[pair][frame] = old[frame][pair]
    for(int p = 0; p < npairs_tot; ++p) {
    	snew(unit_vecs[p], nframes);

    	for(int f = 0; f < nframes; ++f) {
    		copy_rvec(old_unit_vecs[f][p], unit_vecs[p][f]);
    	}
    }

    // Done with old unit vectors matrix
    for(int i = 0; i < nframes; ++i)
        sfree(old_unit_vecs[i]);
    sfree(old_unit_vecs);

    // DEBUG
    // print unit vecs
    FILE *vecf2 = fopen("vecs_pbyf.txt", "w");
    for(int p = 0; p < npairs_tot; ++p) {
        fprintf(vecf2, "\nPair %d, atoms %d and %d:\n", p, corr->atompairs[2 * p], corr->atompairs[2 * p + 1]);
        for(int fr = 0; fr < nframes; ++fr) {
            fprintf(vecf2, "Frame %d: %f, %f, %f\n", fr, unit_vecs[p][fr][XX], unit_vecs[p][fr][YY], unit_vecs[p][fr][ZZ]);
        }
    }
    fclose(vecf2);
    gk_print_log("Unit vectors saved to vecs_pbyf.txt for debugging.\n");


    // calculate autocorrelation from unit vectors

    snew(corr->auto_corr, npairs_tot);

    for(int p = 0; p < npairs_tot; ++p) {
    	snew(corr->auto_corr[p], corr->nt);
    	gc_calc_ac(unit_vecs[p], nframes, corr->nt, corr->auto_corr[p]);
    }

    // TODO: calculate s2

    // Done with unit vectors
    for(int i = 0; i < npairs_tot; ++i)
        sfree(unit_vecs[i]);
    sfree(unit_vecs);
}


void gc_save_corr(struct gcorr_dat_t *corr, const char *corr_fname, const char *s2_fname) {
	if(corr_fname) {
		FILE *f = fopen(corr_fname, "w");
		int p = 0;

		for(int np = 0; np < corr->nnamepairs; ++np) {
			fprintf(f, "# PAIR %s-%s:\n", 
				corr->atomnames[2*np], corr->atomnames[2*np+1]);
			int i;

			for(i = 0; i < corr->natompairs[np]; ++i) {
				fprintf(f, "# Atom %d and %d:\n", 
					corr->atompairs[2*(p+i)], corr->atompairs[2*(p+i)+1]);
				fprintf(f, "# t\tautocorrelation\n");
				fprintf(f, "%f\t%f\n", 
					0.0, 1.0);

				for(int t = 1; t <= corr->nt; ++t) {
					fprintf(f, "%f\t%f\n", 
						t * corr->dt, corr->auto_corr[p+i][t-1]);
				}
			}

			p += i;
		}
		fprintf(f, "\n");
	}

	if(s2_fname) {
		// TODO
	}
}


void gc_free_corr(struct gcorr_dat_t *corr) {
    if(corr->atompairs)   sfree(corr->atompairs);

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
