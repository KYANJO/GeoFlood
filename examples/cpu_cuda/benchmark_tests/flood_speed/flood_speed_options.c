#include "flood_speed_user.h"
#include <fclaw_pointer_map.h>

static int s_user_options_package_id = -1;

static void *
flood_speed_register (user_options_t* user, sc_options_t * opt)
{
    /* [user] User options */
    sc_options_add_bool (opt, 0, "cuda", &user->cuda, 0,
                           "Use cudaclaw [F]");
    sc_options_add_double (opt, 0, "gravity", &user->gravity, 1.0, 
                           "[user] gravity [1.0]");           
    sc_options_add_double (opt, 0, "dry_tolerance", &user->dry_tolerance, 1e-3, 
                           "[user] dry_tolerance [1e-3]");
    sc_options_add_double (opt, 0, "earth_radius", &user->earth_radius, 6367.5e3, 
                           "[user] earth_radius [6367.5e3]");
    sc_options_add_int (opt, 0, "coordinate_system", &user->coordinate_system, 0,
                           "[user] coordinate_system [0]");
    sc_options_add_bool (opt, 0, "mcapa", &user->mcapa, 0,
                           "[user] mcapa [F]");                     

    user->is_registered = 1;
    return NULL;
}

static fclaw_exit_type_t
flood_speed_check (user_options_t *user)
{
    return FCLAW_NOEXIT;
}

static void
flood_speed_destroy(user_options_t *user)
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

    return flood_speed_register(user,opt);
}

static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    user_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (user_options_t*) package;

    return flood_speed_check(user);
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

    flood_speed_destroy (user);

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

user_options_t* flood_speed_options_register (fclaw_app_t * app,
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

user_options_t* flood_speed_get_options(fclaw2d_global_t *glob)
{
    user_options_t* user = (user_options_t*) fclaw_pointer_map_get(glob->options,"flood_speed");
    FCLAW_ASSERT(user != NULL);
    return user;
}

void flood_speed_options_store (fclaw2d_global_t* glob, user_options_t* user)
{
    FCLAW_ASSERT(fclaw_pointer_map_get(glob->options,"flood_speed") == NULL);
    fclaw_pointer_map_insert(glob->options,"flood_speed",user,NULL);
}