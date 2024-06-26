/*
Copyright (c) 2012-2022 Carsten Burstedde, Donna Calhoun, Scott Aiton
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

#ifndef FC2D_CPUCUDA_OPTIONS_H
#define FC2D_CPUCUDA_OPTIONS_H

#include <fclaw_base.h>  /* Needed for fclaw_app_t */

#ifdef __cplusplus
extern "C"
{
#if 0
}
#endif
#endif

struct fclaw2d_global;

typedef struct fc2d_cpucuda_options
{
    // int cuda;
    
    int mwaves;

    const char *order_string;
    int *order;

    int *mthlim;
    const char *mthlim_string;

    int *mthbc;
    const char *mthbc_string;

    int method[7];
    int mcapa;
    int mbathy;
    int src_term;
    int use_fwaves;

    double dry_tolerance_c;
    double wave_tolerance_c;
    int speed_tolerance_entries_c;
    double *speed_tolerance_c;
    const char *speed_tolerance_c_string;

    int ascii_out;  /* Only one type of output now  */    

    int buffer_len;

    int is_registered;
    
} fc2d_cpucuda_options_t;

/**
 * @brief Register options in SC
 * 
 * @param a the app context
 * @param section the section name
 * @param configfile the config file
 * @return fc2d_cpucuda_options_t* a newly allocated options struct
 */
fc2d_cpucuda_options_t *fc2d_cpucuda_options_register (fclaw_app_t * app,
                                                       const char *section,
                                                       const char *configfile);

void fc2d_cpucuda_options_store (struct fclaw2d_global* glob, 
                                 fc2d_cpucuda_options_t* geo_opt);

fc2d_cpucuda_options_t* fc2d_cpucuda_get_options(struct fclaw2d_global *glob);

int cudaclaw_check_parameters(int mwaves);
void cudaclaw_set_method_parameters(int order[], int mthlim[], int mwaves, int use_fwaves);

#ifdef __cplusplus
#if 0
{
#endif
}
#endif

#endif
