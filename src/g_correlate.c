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

#ifdef GTA_BENCH
#include <time.h>
#endif
#include "macros.h"
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

    gmx_bool fft = FALSE;
#ifdef CMEMLIMIT
    gmx_bool mem_limit = TRUE;
#else
    gmx_bool mem_limit = FALSE;
#endif

    init_log("gcorr.log", argc, argv);

    t_filenm fnm[] = {
        {efTRX, "-f", "traj.xtc", ffREAD},
        {efNDX, "-n", "index.ndx", ffOPTRD},
        {efDAT, "-o", "corr.dat", ffWRITE}
    };

    t_pargs pa[] = {
        {"-fft", FALSE, etINT, {&fft}, "use FFT for calculating S2."},
        {"-limit", FALSE, etBOOL, {&mem_limit}, "limit number of trajectory frames loaded into memory. If false, whole trajectory is loaded at once."}
    };

    parse_common_args(&argc, argv, 0, efT_NUMFILES, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

    fnames[efT_TRAJ] = opt2fn("-f", efT_NUMFILES, fnm);
    fnames[efT_NDX] = opt2fn_null("-n", efT_NUMFILES, fnm);
    fnames[efT_OUTDAT] = opt2fn("-o", efT_NUMFILES, fnm);

    // Determine flags
    unsigned long flags = ((int)fft * C_FFT) 
                        | ((int)mem_limit * C_MEM_LIMIT);

    // Run autocorrelation and S2 calculation
    struct corr_dat_t corr;
    calc_ac(fnames, &oenv, &corr, flags);

#ifdef GTA_BENCH
    clock_t clocks = clock() - start;
    print_log("g_tessellate_area took %d clocks, %f seconds.\n", clocks, (float)clocks/CLOCKS_PER_SEC);
#endif

    close_log();

    return 0;
}
