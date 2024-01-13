#ifndef GEOFLOOD_OPTIONS_USER_H
#define GEOFLOOD_OPTIONS_USER_H

#include <fclaw_base.h>  /* Needed for fclaw_app_t */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw2d_global;

typedef struct user_options
{
    int cuda;
    int claw_version;
    int is_registered;

} user_options_t;

/* ------------------------------------- User Options made global ---------------------------------------*/
user_options_t* geoflood_options_register (fclaw_app_t * app,
                                        const char *configfile);

void geoflood_options_store (struct fclaw2d_global *glob, user_options_t* user);

const user_options_t* geoflood_get_options(struct fclaw2d_global *glob);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif