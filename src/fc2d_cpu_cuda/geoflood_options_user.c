#include "geoflood_options_user.h"
#include <fclaw_options.h>
#include <fclaw_pointer_map.h>

#include <fclaw2d_clawpatch_options.h>
#include <fclaw2d_global.h>

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

/* User options made global */
static int s_user_options_package_id = -1;

static void *
geoflood_register (user_options_t *user, sc_options_t * opt)
{
    /* [user] User options */
    sc_options_add_int (opt, 0, "claw-version", &user->claw_version, 5,
                           "Clawpack_version (4 or 5) [5]");

    sc_options_add_bool (opt, 0, "cuda", &user->cuda, 0,
                           "Use cudaclaw5 [F]");

    user->is_registered = 1;

    return NULL;
}

static fclaw_exit_type_t
geoflood_postprocess(user_options_t *user)
{
    /* nothing to post-process yet ... */
    return FCLAW_NOEXIT;
}


static fclaw_exit_type_t
geoflood_check (user_options_t *user)
{
    /* Nothing to check ? */
    return FCLAW_NOEXIT;
}

static void
geoflood_destroy(user_options_t *user)
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

    return geoflood_register(user,opt);
}

static fclaw_exit_type_t
options_postprocess (fclaw_app_t * a, void *package, void *registered)
{
    FCLAW_ASSERT (a != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT (registered == NULL);

    /* errors from the key-value options would have showed up in parsing */
    user_options_t *user = (user_options_t *) package;

    /* post-process this package */
    FCLAW_ASSERT(user->is_registered);

    /* Convert strings to arrays */
    return geoflood_postprocess (user);
}


static fclaw_exit_type_t
options_check(fclaw_app_t *app, void *package,void *registered)
{
    user_options_t           *user;

    FCLAW_ASSERT (app != NULL);
    FCLAW_ASSERT (package != NULL);
    FCLAW_ASSERT(registered == NULL);

    user = (user_options_t*) package;

    return geoflood_check(user);
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

    geoflood_destroy (user);

    FCLAW_FREE (user);
}


static const fclaw_app_options_vtable_t options_vtable_user =
{
    options_register,
    options_postprocess,
    options_check,
    options_destroy
};

/* --------------------- Public interface access functions ---------------------------- */

user_options_t* geoflood_options_register (fclaw_app_t * app, 
                                        const char *section,
                                        const char *configfile)
{
    user_options_t *user;
    FCLAW_ASSERT (app != NULL);

    user = FCLAW_ALLOC (user_options_t, 1);
    fclaw_app_options_register (app, section, configfile, 
                                &options_vtable_user, user);

    fclaw_app_set_attribute(app,section,user);
    return user;
}

user_options_t* geoflood_get_options(fclaw2d_global_t* glob)
{
    // int id = s_user_options_package_id;
    // return (user_options_t*) fclaw_package_get_options(glob, id);    
    user_options_t* user = (user_options_t*) 
	   							fclaw_pointer_map_get(glob->options, "fc2d_cpu_cuda");
	FCLAW_ASSERT(user != NULL);
	return user;
}

void geoflood_options_store (fclaw2d_global_t* glob, user_options_t* user)
{
    // FCLAW_ASSERT(s_user_options_package_id == -1);
    // int id = fclaw_package_container_add_pkg(glob,user);
    // s_user_options_package_id = id;
    FCLAW_ASSERT(fclaw_pointer_map_get(glob->options,"fc2d_cpu_cuda") == NULL);
	fclaw_pointer_map_insert(glob->options, "fc2d_cpu_cuda", user, NULL);
}

#ifdef __cplusplus
#if 0
{
#endif
}
#endif