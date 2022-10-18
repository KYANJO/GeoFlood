/*
# ----------------------------------------------
# @author:  Brian Kyanjo
# @contact: briankyanjo@u.boisestate.edu
# @date:    2022-10-16
# @version: 1.0
# @desc:    
# ------------------------------------------------
*/

#include "malpasset_user.h"
#include <fclaw2d_include_all.h>

#include <fclaw2d_clawpatch.h>
#include <fclaw2d_clawpatch_options.h>

#include <fc2d_geoclaw.h>
#include <fc2d_geoclaw_options.h>


static
fclaw2d_domain_t* create_domain(sc_MPI_Comm mpicomm, 
                                fclaw_options_t* gparms)
{
    p4est_connectivity_t     *conn = NULL;
    fclaw2d_domain_t         *domain;
    fclaw2d_map_context_t    *cont = NULL;

    /* Size is set by [ax,bx] x [ay, by], set in .ini file */
    conn = p4est_connectivity_new_unitsquare();
    cont = fclaw2d_map_new_nomap();

    domain = fclaw2d_domain_new_conn_map (mpicomm, gparms->minlevel, conn, cont);
    fclaw2d_domain_list_levels(domain, FCLAW_VERBOSITY_ESSENTIAL);
    fclaw2d_domain_list_neighbors(domain, FCLAW_VERBOSITY_DEBUG);

    return domain;
}

static
void run_program(fclaw2d_global_t* glob)
{
    fclaw2d_domain_t    **domain = &glob->domain;

    fclaw2d_domain_data_new(*domain);

    fclaw2d_vtables_initialize(glob);

    fc2d_geoclaw_solver_initialize();

    /* ---------------------------------------------------------------
       Run
       --------------------------------------------------------------- */
    fc2d_geoclaw_module_setup(glob);

    fclaw2d_initialize(glob);
    fc2d_geoclaw_run(glob);
    
    fclaw2d_finalize(glob);
}

int
main (int argc, char **argv)
{
    fclaw_app_t *app;
    int first_arg;
    fclaw_exit_type_t vexit;

    /* Options */
    sc_options_t                *options;
    fclaw_options_t             *gparms;
    fclaw2d_clawpatch_options_t *clawpatchopt;
    fc2d_geoclaw_options_t      *geoclawopt;

    sc_MPI_Comm mpicomm;
    fclaw2d_domain_t* domain;
    fclaw2d_global_t* glob;

    int retval;

     /* Initialize application */
    app = fclaw_app_new (&argc, &argv, NULL);

    gparms                   = fclaw_options_register(app, NULL,"fclaw_options.ini");
    clawpatchopt = fclaw2d_clawpatch_options_register(app,"clawpatch", "fclaw_options.ini");
    geoclawopt        = fc2d_geoclaw_options_register(app, "geoclaw","fclaw_options.ini");

     /* Read configuration file(s) and command line, and process options */
    options = fclaw_app_get_options (app);
    retval = fclaw_options_read_from_file(options);
    vexit =  fclaw_app_options_parse (app, &first_arg,"fclaw_options.ini.used");

    if (!retval & !vexit)
    {
        mpicomm = fclaw_app_get_mpi_size_rank (app, NULL, NULL);
        domain = create_domain(mpicomm, gparms);
    
        glob = fclaw2d_global_new();
        fclaw2d_global_store_domain(glob, domain);

        fclaw2d_options_store           (glob, gparms);
        fclaw2d_clawpatch_options_store (glob, clawpatchopt);
        fc2d_geoclaw_options_store      (glob, geoclawopt);

        /* Run the program */
        run_program(glob);

        fclaw2d_global_destroy(glob);
    }

    fclaw_app_destroy (app);

    return 0;

}
