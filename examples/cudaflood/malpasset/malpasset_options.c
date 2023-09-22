/*
Copyright (c) 2012-2021 Carsten Burstedde, Donna Calhoun
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

 * Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "malpasset_user.h"

static int s_user_options_package_id = -1;

static void *
malpasset_register (user_options_t* user, sc_options_t * opt)
{
    /* [user] User options */
    sc_options_add_bool (opt, 0, "cuda", &user->cuda, 0,
                           "Use cudaclaw [F]");    
    sc_options_add_double (opt, 0, "gravity", &user->gravity, 1.0, 
                           "[user] gravity [1.0]");
    sc_options_add_double (opt, 0, "dry_tolerance", &user->dry_tolerance, 0.001, 
                           "[user] dry_tolerance [0.001]");
    sc_options_add_double (opt, 0, "earth_radius", &user->earth_radius, 6367.5e3, 
                           "[user] earth_radius [6367.5e3]");  
    sc_options_add_int (opt, 0, "coordinate_system", &user->coordinate_system, 0,
                           "[user] coordinate_system [0]");
    sc_options_add_int (opt, 0, "mcapa", &user->mcapa, 0,
                           "[user] mcapa [0]");                            

    user->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
malpasset_check (user_options_t *user)
{
    if(user->cuda != 0)
    {
        if (user->example != 0) {
            fclaw_global_essentialf ("Option --user:example must be 0\n");
            return FCLAW_EXIT_QUIET;
        }
    }
    return FCLAW_NOEXIT;
}

static void
malpasset_destroy(user_options_t *user)
{
    /* Nothing to destroy */
}

/* ------- Generic option handling routines that call above routines ----- */
static void*
options_register (fclaw_app_t * app, void *package, sc_options_t * opt)
{
    user_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (opt != NULL);

    user = (user_options_t*) package;

    return malpasset_register(user,opt);
}

static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    user_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (user_options_t*) package;

    return malpasset_check(user);
}


static void
options_destroy (fclaw_app_t * app, void *package, void *registered)
{
    user_options_t *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    user = (user_options_t*) package;
    FCLAW_ASSERT (user->is_registered);

    malpasset_destroy (user);

    FCLAW_FREE (user);
}


static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    NULL,
    options_check,
    options_destroy
};

/* ------------- User options access functions --------------------- */

user_options_t* malpasset_options_register (fclaw_app_t * app,
                                          const char *configfile)
{
    user_options_t *user;
    FCLAW_ASSERT (app != NULL);

    user = FCLAW_ALLOC (user_options_t, 1);
    fclaw_app_options_register (app,"user", configfile, &options_vtable_user,
                                user);

    fclaw_app_set_attribute(app,"user",user);
    return user;
}

void malpasset_options_store (fclaw2d_global_t* glob, user_options_t* user)
{
    FCLAW_ASSERT(s_user_options_package_id == -1);
    int id = fclaw_package_container_add_pkg(glob,user);
    s_user_options_package_id = id;
}

user_options_t* malpasset_get_options(fclaw2d_global_t* glob)
{
    int id = s_user_options_package_id;
    return (user_options_t*) fclaw_package_get_options(glob, id);    
}
