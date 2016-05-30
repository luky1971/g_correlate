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

#include <string.h> // for strlen
#ifdef GTA_BENCH
#include <time.h>
#endif
#include "macros.h" // for asize(), used with parse_common_args()
#include "smalloc.h" // memory allocation
#include "gkut_io.h" // for topology i/o
#include "gkut_log.h"
#include "correlate.h"


int main(int argc, char *argv[]) {
#ifdef GTA_BENCH
    clock_t start = clock();
#endif

    const char *desc[] = {
        "g_correlate calculates autocorrelation functions and generalized (S2) order parameters for MD trajectories. \n",
        "It reads in a trajectory file through the -f option (supported formats=xtc,trr,pdb). \n"
    };

    const char *fnames[efT_NUMFILES];
    output_env_t oenv = NULL;

    struct gcorr_dat_t corr;
    gc_init_corr_dat(&corr); // set default values

    char *pairnames = NULL;
    gmx_bool fft = FALSE;
    gmx_bool mem_limit;
#ifdef CMEMLIMIT
    mem_limit = TRUE;
#else
    mem_limit = FALSE;
#endif

    // Initialize log file
    gk_init_log("gcorr.log", argc, argv);

    t_filenm fnm[] = {
        {efTRX, "-f", "traj.xtc", ffREAD},
        {efNDX, "-n", "index.ndx", ffOPTRD},
        {efSTX, "-s", "topol.tpr", ffREAD},
        {efDAT, "-o", "corr.dat", ffWRITE},
        {efDAT, "-s2", "s2.dat", ffWRITE}
    };

    t_pargs pa[] = {
        {"-a", FALSE, etSTR, {&pairnames}, 
            "Comma-delimited list of atom name pairs to use for calculating S2, supports wildcards (ex. 'N-H,ND2-H*,C*-H*')"},
        {"-dt", FALSE, etREAL, {&(corr.dt)}, 
            "The time delay step for the domain of the autocorrelation functions. The default is to use the trajectory's timestep."},
        {"-nt", FALSE, etINT, {&(corr.nt)}, 
            "The number of time delays in the domain of the autocorrelation functions. The default is for nt*dt to go up to the length of the trajectory."},
        // {"-fft", FALSE, etINT, {&fft}, "Use FFT for calculating S2."},
        // {"-limit", FALSE, etBOOL, {&mem_limit}, 
                // "Limit number of trajectory frames loaded into memory. If false, whole trajectory is loaded at once."}
    };

    parse_common_args(&argc, argv, 0, efT_NUMFILES, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

    if(pairnames == NULL)
        gk_log_fatal(FARGS, "Must provide atom name pairs using -a option!\n");

    // Get filenames
    fnames[efT_TRAJ] = opt2fn("-f", efT_NUMFILES, fnm);
    fnames[efT_NDX] = opt2fn_null("-n", efT_NUMFILES, fnm);
    fnames[efT_TOP] = opt2fn("-s", efT_NUMFILES, fnm);
    fnames[efT_OUTDAT] = opt2fn("-o", efT_NUMFILES, fnm);
    fnames[efT_S2DAT] = opt2fn("-s2", efT_NUMFILES, fnm);

    // Set flags
    unsigned long flags = ((int)fft * C_FFT) 
                        | ((int)mem_limit * C_MEM_LIMIT);


    // Initialize autocorrelation parameters
    snew(corr.atomnames, strlen(pairnames));

    // Parse atom name pairs for autocorrelator input
    int nwords = 0;
    inparse(pairnames, corr.atomnames, &nwords);

    if(nwords <= 0)
        gk_log_fatal(FARGS, "No atom pairs specified using -a option!\n");
    if(nwords % 2 != 0)
        gk_log_fatal(FARGS, "A specified atom name (in option -a) is unpaired!\n");

    srenew(corr.atomnames, nwords); // Free excess memory
    corr.nnamepairs = nwords / 2;

    // DEBUG
    // printf("%d pairs parsed.\n", corr.nnamepairs);
    // for(int i = 0; i < corr.nnamepairs; ++i) {
    //     printf("%s and %s\n", corr.atomnames[2*i], corr.atomnames[2*i+1]);
    // }


    // Run autocorrelation and S2 calculation
    gc_correlate(fnames, &oenv, &corr, flags);


    // Print results
    t_topology top;
    t_atoms *pAtoms = NULL;

    if(gk_read_topology(fnames[efT_TOP], &top))
        pAtoms = &(top.atoms);

    gc_save_corr(&corr, fnames[efT_OUTDAT], pAtoms);
    gc_save_s2(&corr, fnames[efT_S2DAT], pAtoms);

    gk_print_log("Bye!\n");


    // Cleanup
    if(pAtoms)
        gk_free_topology(&top);
    gc_free_corr(&corr);
    sfree(corr.atomnames);

#ifdef GTA_BENCH
    clock_t clocks = clock() - start;
    gk_print_log("g_tessellate_area took %d clocks, %f seconds.\n", clocks, (float)clocks/CLOCKS_PER_SEC);
#endif

    gk_close_log();

    return 0;
}
