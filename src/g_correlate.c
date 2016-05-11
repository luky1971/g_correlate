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

#include <ctype.h> // for isspace
#include <string.h> // for strlen
#ifdef GTA_BENCH
#include <time.h>
#endif
#include "macros.h" // for asize(), used with parse_common_args()
#include "smalloc.h" // memory allocation
#include "gkut_log.h"
#include "correlate.h"

// NOTE: this modifies the inwords array and simply sets pointers to different portions of it in outwords,
// to avoid allocating and copying memory.
static void inparse(char *inwords, const char **outwords, int *nwords) {
    char *c = inwords;
    const char **w = outwords;
    int newword = 1;
    while(c != NULL && *c != '\0') {
        if(isspace(*c) || *c == '-' || *c == ',') {
            *c = '\0';
            newword = 1;
        }
        else if(newword) {
            *(w++) = c;
            newword = 0;
        }
        ++c;
    }
    *nwords = w - outwords;
}

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
        {efDAT, "-o", "corr.dat", ffWRITE}
    };

    t_pargs pa[] = {
        {"-a", FALSE, etSTR, {&pairnames}, "Comma-delimited list of atom name pairs to use for calculating S2, supports wildcards (ex. 'N-H,ND2-H*,C*-H*')"},
        {"-fft", FALSE, etINT, {&fft}, "Use FFT for calculating S2."},
        {"-limit", FALSE, etBOOL, {&mem_limit}, "Limit number of trajectory frames loaded into memory. If false, whole trajectory is loaded at once."}
    };

    parse_common_args(&argc, argv, 0, efT_NUMFILES, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv);

    if(pairnames == NULL)
        gk_log_fatal(FARGS, "Must provide atom name pairs using -a option!\n");

    // Get filenames
    fnames[efT_TRAJ] = opt2fn("-f", efT_NUMFILES, fnm);
    fnames[efT_NDX] = opt2fn_null("-n", efT_NUMFILES, fnm);
    fnames[efT_TOP] = opt2fn("-s", efT_NUMFILES, fnm);
    fnames[efT_OUTDAT] = opt2fn("-o", efT_NUMFILES, fnm);

    // Set flags
    unsigned long flags = ((int)fft * C_FFT) 
                        | ((int)mem_limit * C_MEM_LIMIT);


    // Initialize autocorrelation parameters
    
    struct corr_dat_t corr;
    snew(corr.atomnames, strlen(pairnames));

    // Parse atom name pairs for autocorrelator input
    int nwords = 0;
    inparse(pairnames, corr.atomnames, &nwords);

    if(nwords <= 0)
        gk_log_fatal(FARGS, "No atom pairs specified using -a option!\n");
    if(nwords % 2 != 0)
        gk_log_fatal(FARGS, "A specified atom name (in option -a) is unpaired!\n");

    srenew(corr.atomnames, nwords); // Free excess memory
    corr.npairs = nwords / 2;

    // DEBUG
    // printf("%d pairs parsed.\n", corr.npairs);
    // for(int i = 0; i < corr.npairs; ++i) {
    //     printf("%s and %s\n", corr.atomnames[2*i], corr.atomnames[2*i+1]);
    // }

    // Run autocorrelation and S2 calculation
    calc_ac(fnames, &oenv, &corr, flags);

    // TODO: calculate S2, print results, and whatnot!

    // Cleanup
    free_corr(&corr);
    sfree(corr.atomnames);

#ifdef GTA_BENCH
    clock_t clocks = clock() - start;
    gk_print_log("g_tessellate_area took %d clocks, %f seconds.\n", clocks, (float)clocks/CLOCKS_PER_SEC);
#endif

    gk_close_log();

    return 0;
}
