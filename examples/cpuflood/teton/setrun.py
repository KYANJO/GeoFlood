# ----------------------------------------------
# @author:  Brian Kyanjo
# @contact: briankyanjo@u.boisestate.edu
# @date:    2022-10-16
# @version: 1.0
# ------------------------------------------------

import os
import sys
import numpy as np
from pdb import *

import tools

#===============================================================================
# Importing scripts dictionary
#===============================================================================
sys.path.append('../../../scripts')
import geoflood # -- importing geoflood.py

#===============================================================================
# scratch directory
#===============================================================================
scratch_dir = os.path.join('../scratch')

#===============================================================================
# User specified parameters
#===============================================================================
#------------------ Time stepping------------------------------------------------
initial_dt = 0.0001  # Initial time step
fixed_dt = False   # Take constant time step

# -------------------- Output files -------------------------------------------------
output_style = 1 #changed 10.21   

if output_style == 1:
    # Total number of frames will be frames_per_minute*60*n_hours

    n_hours = 15              # Total number of hours in simulation, changed 10.14.2020  should be 5      
    

    frames_per_minute = 1/30   # Frames every 1/2 hour

if output_style == 2:
    output_times = [1,2,3]    # Specify exact times to output files

if output_style == 3:
    step_interval = 10   # Create output file every 10 steps
    total_steps = 500    # ... for a total of 500 steps (so 50 output files total)

#-------------------  Computational coarse grid ---------------------------------------
mx =  32
my = 32

minlevel = 2
maxlevel = 5 #resolution based on levels
ratios_x = [2,4,4,4]
ratios_y = [2,4,4,4]
ratios_t = [2,4,4,4] #should this be 0,0,0,0?

#-------------------manning coefficient -----------------------------------------------
manning_coefficient = 0.06

#-------------------  Number of dimensions ---------------------------------------
num_dim = 2

# --------------------- Topography file -----------------------------------------------
# topofile = 'topos/TetonLarge.topo'
topofile = 'topos/TetonDamLatLong.topo'

#------------------------------
def setrun(claw_pkg='geoclaw'):
#------------------------------

    """
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "geoclaw" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData

    """

    from clawpack.clawutil import data

    assert claw_pkg.lower() == 'geoclaw',  "Expected claw_pkg = 'geoclaw'"

   

    rundata = data.ClawRunData(claw_pkg, num_dim)

    #------------------------------------------------------------------
    # GeoClaw specific parameters:
    #------------------------------------------------------------------
    rundata = setgeo(rundata)

    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #   (or to amr2ez.data for AMR)
    #------------------------------------------------------------------
    clawdata = rundata.clawdata  # initialized when rundata instantiated


    # Set single grid parameters first.
    # See below for AMR parameters.


    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.num_dim = num_dim

    m_topo,n_topo,xllcorner,yllcorner,cellsize = tools.read_topo_data(topofile)


    # Topo info (TetonDamLatLong.topo)
    # m_topo = 4180
    # n_topo = 1464
    # xllcorner = -112.390734400000
    # yllcorner = 43.581746970335
    # cellsize = 0.000277729665

    # Derived info from the topo map
    mx_topo = m_topo - 1
    my_topo = n_topo - 1
    xurcorner = xllcorner + cellsize*mx_topo
    yurcorner = yllcorner + cellsize*my_topo
    ll_topo = np.array([xllcorner, yllcorner])
    ur_topo = np.array([xurcorner, yurcorner])

    print("")
    print("Topo domain")
    print("%-12s (%14.8f, %12.8f)" % ("Lower left",ll_topo[0],ll_topo[1]))
    print("%-12s (%14.8f, %12.8f)" % ("Upper right",ur_topo[0],ur_topo[1]))
    print("")

    dims_topo = ur_topo - ll_topo

    # Try to match aspect ratio of topo map
    clawdata.num_cells[0] = mx
    clawdata.num_cells[1] = my

    dim_topo = ur_topo - ll_topo
    mdpt_topo = ll_topo + 0.5*dim_topo

    dim_comp = 0.975*dim_topo   # Shrink domain inside of given bathymetry.

    clawdata.lower[0] = mdpt_topo[0] - dim_comp[0]/2.0
    clawdata.upper[0] = mdpt_topo[0] + dim_comp[0]/2.0

    clawdata.lower[1] = mdpt_topo[1] - dim_comp[1]/2.0
    clawdata.upper[1] = mdpt_topo[1] + dim_comp[1]/2.0

    print("")
    print("Computational domain")
    print("%-12s (%14.8f, %12.8f)" % ("Lower left",clawdata.lower[0],clawdata.lower[1]))
    print("%-12s (%14.8f, %12.8f)" % ("Upper right",clawdata.upper[0],clawdata.upper[1]))
    print("")

    print("Approximate aspect ratio : {0:16.8f}".format(float(clawdata.num_cells[0])/clawdata.num_cells[1]))
    print("Actual      aspect ratio : {0:16.8f}".format(dims_topo[0]/dims_topo[1]))

    # print "[{0:20.12f},{1:20.12f}]".format(*clawdata.lower)
    # print "[{0:20.12f},{1:20.12f}]".format(*clawdata.upper)

    dims_computed = np.array([clawdata.upper[0]-clawdata.lower[0], clawdata.upper[1]-clawdata.lower[1]])
    print("Computed aspect ratio    : {0:20.12f}".format(dims_computed[0]/dims_computed[1]))

    print("")
    print("Details in km : ")    

    lon = np.array([clawdata.lower[0],clawdata.upper[0]])
    lat = np.array([clawdata.lower[1],clawdata.upper[1]])
    d = tools.compute_distances(lon,lat)
   
    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.num_eqn = 3

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.num_aux = 3

    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.capa_index = 2

    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = 0.0


    # Restart from checkpoint file of a previous run?
    # Note: If restarting, you must also change the Makefile to set:
    #    RESTART = True
    # If restarting, t0 above should be from original run, and the
    # restart_file 'fort.chkNNNNN' specified below should be in
    # the OUTDIR indicated in Makefile.

    clawdata.restart = False               # True to restart from prior results
    clawdata.restart_file = 'fort.chk00006'  # File to use for restart data

    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.output_style = output_style

    if clawdata.output_style == 1:
        # Output nout frames at equally spaced times up to tfinal:
        # n_hours = 20.0
        # frames_per_minute = 1/30.0 # Frames every 5 seconds
        clawdata.num_output_times = int(frames_per_minute*60*n_hours)  # Plot every 10 seconds
        clawdata.tfinal = 60*60*n_hours
        clawdata.output_t0 = True  # output at initial (or restart) time?

    elif clawdata.output_style == 2:
        # Specify a list of output times.
        clawdata.output_times = [0.5, 1.0]

    elif clawdata.output_style == 3:
        # Output every iout timesteps with a total of ntot time steps:
        clawdata.output_step_interval = step_interval
        clawdata.total_steps = total_steps
        clawdata.output_t0 = True
        clawdata.tfinal = total_steps*fixed_dt

    clawdata.output_format = 'ascii'      # 'ascii' or 'netcdf'

    clawdata.output_q_components = 'all'   # could be list such as [True,True]
    clawdata.output_aux_components = 'none'  # could be list
    clawdata.output_aux_onlyonce = True    # output aux arrays only at t0



    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = 1



    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = not fixed_dt

    # Initial time step for variable dt.
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = initial_dt

    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = 1e+99

    # Desired Courant number if variable dt used, and max to allow without
    # retaking step with a smaller dt:
    clawdata.cfl_desired = 0.75
    clawdata.cfl_max = 1.0

    # Maximum number of time steps to allow between output times:
    clawdata.steps_max = 5000

    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = 2

    # Use dimensional splitting? (not yet available for AMR)
    clawdata.dimensional_split = 'unsplit'

    # For unsplit method, transverse_waves can be
    #  0 or 'none'      ==> donor cell (only normal solver used)
    #  1 or 'increment' ==> corner transport of waves
    #  2 or 'all'       ==> corner transport of 2nd order corrections too
    clawdata.transverse_waves = 2

    # Number of waves in the Riemann solution:
    clawdata.num_waves = 3

    # List of limiters to use for each wave family:
    # Required:  len(limiter) == num_waves
    # Some options:
    #   0 or 'none'     ==> no limiter (Lax-Wendroff)
    #   1 or 'minmod'   ==> minmod
    #   2 or 'superbee' ==> superbee
    #   3 or 'mc'       ==> MC limiter
    #   4 or 'vanleer'  ==> van Leer
    clawdata.limiter = ['mc', 'mc', 'mc']

    clawdata.use_fwaves = True    # True ==> use f-wave version of algorithms

    # Source terms splitting:
    #   src_split == 0 or 'none'    ==> no source term (src routine never called)
    #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used,
    #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.
    clawdata.source_split = 'godunov'


    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.num_ghost = 2

    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity

    clawdata.bc_lower[0] = 'extrap'
    clawdata.bc_upper[0] = 'extrap'

    clawdata.bc_lower[1] = 'extrap'
    clawdata.bc_upper[1] = 'extrap'

    # Specify when checkpoint files should be created that can be
    # used to restart a computation.

    clawdata.checkpt_style = 0

    if clawdata.checkpt_style == 0:
        # Do not checkpoint at all
        pass

    elif np.abs(clawdata.checkpt_style) == 1:
        # Checkpoint only at tfinal.
        pass

    elif np.abs(clawdata.checkpt_style) == 2:
        # Specify a list of checkpoint times.  
        clawdata.checkpt_times = [0.1,0.15]

    elif np.abs(clawdata.checkpt_style) == 3:
        # Checkpoint every checkpt_interval timesteps (on Level 1)
        # and at the final time.
        clawdata.checkpt_interval = 5

    # --------------------------------------------------------
    # GeoFlood parameters. 
    # These will overwrite any similar parameters listed above
    # --------------------------------------------------------

    geoflooddata = geoflood.GeoFlooddata()
    geoflooddata.minlevel = minlevel
    geoflooddata.maxlevel = maxlevel

    geoflooddata.regrid_interval = 3
    geoflooddata.refine_threshold = 0.01
    geoflooddata.coarsen_threshold = 0.005

    geoflooddata.subcycle = False
    geoflooddata.output = True
    geoflooddata.output_gauges = True

    # geoflood verbosity choices : 
    # 0 or 'silent'      : No output to the terminal
    # 1 or 'essential'   : Only essential output, including errors.
    # 2 or 'production'  : Production level output
    # 3 or 'info'        : More detailed output
    # 4 or 'debug'       : Includes detailed output from each processor
    geoflooddata.verbosity = 'production'


    # Block dimensions for non-square domains
    geoflooddata.mi = 5
    geoflooddata.mj = 2

    geoflooddata.user = {'example'     : 1}

    # -----------------------------------------------
    # AMR parameters:
    # -----------------------------------------------
    amrdata = rundata.amrdata

    amrdata.amr_levels_max = maxlevel    # Set to 3 for best results

    amrdata.refinement_ratios_x = ratios_x 
    amrdata.refinement_ratios_y = ratios_y
    amrdata.refinement_ratios_t = ratios_t

    # rundata.tol = -1
    # rundata.tolsp = 0.001

    # Specify type of each aux variable in amrdata.auxtype.
    # This must be a list of length maux, each element of which is one of:
    #   'center',  'capacity', 'xleft', or 'yleft'  (see documentation).

    amrdata.aux_type = ['center','capacity','yleft','center']


    # Flag using refinement routine flag2refine rather than richardson error
    amrdata.flag_richardson = False    # use Richardson?
    amrdata.flag2refine = True
    amrdata.regrid_interval = 3
    amrdata.regrid_buffer_width  = 2
    amrdata.clustering_cutoff = 0.700000
    amrdata.verbosity_regrid = 0

    # -----------------------------------------------
    # INL Regions
    #   Regions to be refined :
    #    (1) Refine initial reservoir to level 4
    #        (otherwise, we won't resolve valley, and
    #        won't fill the reservoir properly)
    #    (2) Refine around nuclear power plant (indicated by gauge
    #        100, 101, ..., 115, below)
    #    (3) Computational domain, with maxlevel=4
    #
    # To specify regions of refinement append lines of the form
    #    regions.append([minlevel,maxlevel,t1,t2,x1,x2,y1,y2])

    # -----------------------------------------------
    regions = rundata.regiondata.regions

    # Region containing initial reservoir
    regions.append([maxlevel,maxlevel, 0, 1.e10,-111.543,-111.24,43.88, 43.965])

    # Box containing gauge location locations
    xll = [-111.64, 43.913661]  # From email
    xur = [-111.60, 43.92]  # from email
    region_lower, region_upper,_ = tools.region_coords(xll,xur,
                                                    clawdata.num_cells,
                                                    clawdata.lower,
                                                    clawdata.upper)

    regions.append([maxlevel-1,maxlevel-1,0, 1e10, region_lower[0],region_upper[0],
                    region_lower[1],region_upper[1]])

    # Computational domain.  With exception of region above, don't go beyond level 4
    regions.append([minlevel-1,minlevel+1,0, 1e10, clawdata.lower[0],clawdata.upper[0],
                    clawdata.lower[1],clawdata.upper[1]])

    # -------------------------------------------------------
    # INL Gauges
    #     -- Set gauges at Teton City and Wilford
    #     -- Remaining gauges build border around power plant
    #
    # For gauges append lines of the form  [gaugeno, x, y, t1, t2]
    # -------------------------------------------------------
    rundata.gaugedata.gtype = {}
   
    #Stationary Gauges

    #Teton_Canyon_
    xc,yc = [-111.593965, 43.934059] 
    rundata.gaugedata.gauges.append([1,xc,yc,0.,clawdata.tfinal])  # Mid Teton Canyon 
    rundata.gaugedata.gtype[1] = 'stationary'

    #Teton_Canyon_Mouth_
    xc,yc = [-111.66637, 43.933847] 
    rundata.gaugedata.gauges.append([2,xc,yc,0.,clawdata.tfinal])  # Teton Canyon Mouth 
    rundata.gaugedata.gtype[2] = 'stationary'

    #Wilford_Gauge_
    xc,yc = [-111.672, 43.9144]
    rundata.gaugedata.gauges.append([3,xc,yc,0.,clawdata.tfinal])  # Wilford Gauge 
    rundata.gaugedata.gtype[3] = 'stationary'

    #Sugar_City_Gauge_
    xc,yc = [-111.743358, 43.873840]
    rundata.gaugedata.gauges.append([4,xc,yc,0.,clawdata.tfinal])  # Sugar City Gauge 2 
    rundata.gaugedata.gtype[4] = 'stationary'

    #Roberts Gauge 
    xc,yc = [-112.126403, 43.7202] 
    rundata.gaugedata.gauges.append([5,xc,yc,0.,clawdata.tfinal])  # Roberts Gauge    
    rundata.gaugedata.gtype[5] = 'stationary'

    #Rexburg_Gauge_
    xc,yc = [-111.792295, 43.823048] 
    rundata.gaugedata.gauges.append([6,xc,yc,0.,clawdata.tfinal])  # Rexburg Gauge 
    rundata.gaugedata.gtype[6] = 'stationary'

    # or to have some of each type, use a dictionary:
    rundata.gaugedata.gtype = {}
    
    # lagrangian gauges Northeastern
    for iyg in range(0,3): #ten is the grid
        for ixg in range(0,3):
            gaugeno = 10*iyg + ixg + 100
            yg = 43.91335 + 0.01*iyg #testing
            xg = -111.6211080 + 0.012*ixg #test
            rundata.gaugedata.gauges.append([gaugeno, xg, yg, 0., 1e10])
            rundata.gaugedata.gtype[gaugeno] = 'lagrangian'
    
    # lagrangian gauges Southwestern - Menan Butte
    for iyg in range(0,3): #ten is the grid
        for ixg in range(0,3):
            gaugeno = 10*iyg + ixg + 200
            yg = 43.800575 + 0.01*iyg #testing
            xg = -111.9419740 + 0.012*ixg #test
            rundata.gaugedata.gauges.append([gaugeno, xg, yg, 0., 1e10])
            rundata.gaugedata.gtype[gaugeno] = 'lagrangian'
    # # Wilford
    # xc,yc = [-111.672222,43.914444]
    # rundata.gaugedata.gauges.append([1,xc,yc,0.,clawdata.tfinal])  # Wilford

    # # Teton City
    # xc,yc = [-111.669167,43.887778]
    # rundata.gaugedata.gauges.append([2,xc,yc,0.,clawdata.tfinal])  # Teton City

    # Power plant, with border constructed of 4*m gauges
    # Start at SW corner; build gauges in counter-clockwise order in a
    # square around the region [xll,xur].
    #-------------------------------------------------------
    # m = 2  # Gauge spacing along one edge (m=4 --> edge divided into four sections)
    # gauge_counter = 100

    # # South West corner of power plant
    # xll = [-111.623926, 43.913661]  # From email

    # # North East corner of power plant
    # xur = [-111.620150, 43.916382]  # from email

    # s = np.linspace(0,1.,m+1)
    # for i in range(0,m):
    #     x = xll[0] + (xur[0] - xll[0])*s[i]
    #     rundata.gaugedata.gauges.append([gauge_counter,x,xll[1],0.,clawdata.tfinal])
    #     gauge_counter = gauge_counter + 1

    # for i in range(0,m):
    #     y = xll[1] + (xur[1] - xll[1])*s[i]
    #     rundata.gaugedata.gauges.append([gauge_counter,xur[0],y,0.,clawdata.tfinal])
    #     gauge_counter = gauge_counter + 1

    # for i in range(0,m):
    #     x = xur[0] + (xll[0] - xur[0])*s[i]
    #     rundata.gaugedata.gauges.append([gauge_counter,x,xur[1],0.,clawdata.tfinal])
    #     gauge_counter = gauge_counter + 1

    # for i in range(0,m):
    #     y = xur[1] + (xll[1] - xur[1])*s[i]
    #     rundata.gaugedata.gauges.append([gauge_counter,xll[0],y,0.,clawdata.tfinal])
    #     gauge_counter = gauge_counter + 1


    # -------------------------------------------------------
    # For developers
    #    -- Toggle debugging print statements:
    # -------------------------------------------------------
    amrdata.dprint = False      # print domain flags
    amrdata.eprint = False      # print err est flags
    amrdata.edebug = False      # even more err est flags
    amrdata.gprint = False      # grid bisection/clustering
    amrdata.nprint = False      # proper nesting output
    amrdata.pprint = False      # proj. of tagged points
    amrdata.rprint = False      # print regridding summary
    amrdata.sprint = False      # space/memory output
    amrdata.tprint = True      # time step reporting each level
    amrdata.uprint = False      # update/upbnd reporting


    return rundata,geoflooddata
    # end of function setrun
    # ----------------------

#-------------------
def setgeo(rundata):
#-------------------
    """
    Set GeoClaw specific runtime parameters.
    For documentation see ....
    """

    try:
        geo_data = rundata.geo_data
    except:
        print("*** Error, this rundata has no geo_data attribute")
        raise AttributeError("Missing geo_data attribute")

    # == Physics ==
    geo_data.gravity = 9.81
    geo_data.coordinate_system = 2   # LatLong coordinates
    geo_data.earth_radius = 6367.5e3

    # == Forcing Options
    geo_data.coriolis_forcing = True

    # == Algorithm and Initial Conditions ==
    geo_data.sea_level = 0.0
    geo_data.dry_tolerance = 1.e-3
    geo_data.friction_forcing = True
    geo_data.manning_coefficient = manning_coefficient
    geo_data.friction_depth = 1.e6

    # Refinement data
    refinement_data = rundata.refinement_data
    refinement_data.wave_tolerance = 1.e-2
    refinement_data.deep_depth = 1e2
    refinement_data.max_level_deep = 3
    refinement_data.variable_dt_refinement_ratios = False

    # == settopo.data values ==
    topo_data = rundata.topo_data
    # for topography, append lines of the form
    #    [topotype, minlevel, maxlevel, t1, t2, fname]

    topo_data.topofiles.append([2, 1, 10, 0, 1e10, topofile])


    # == setdtopo.data values ==
    topo_data = rundata.topo_data
    # for moving topography, append lines of the form :   (<= 1 allowed for now!)
    #   [topotype, minlevel,maxlevel,fname]

    # == setqinit.data values ==
    rundata.qinit_data.qinit_type = 0
    rundata.qinit_data.qinitfiles = []
    # for qinit perturbations, append lines of the form: (<= 1 allowed for now!)
    #   [minlev, maxlev, fname]

    # == setfixedgrids.data values ==
    # fixedgrids = rundata.fixed_grid_data
    # fixedgrids = rundata.fixed_grid_data.fixedgrids
    # for fixed grids append lines of the form
    # [t1,t2,noutput,x1,x2,y1,y2,xpoints,ypoints,\
    #  ioutarrivaltimes,ioutsurfacemax]

    return rundata
    # end of function setgeo
    # ----------------------

if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    rundata,geoflooddata = setrun(*sys.argv[1:])
    rundata.write()

    geoflooddata.write(rundata)  # writes a geoflood geoflood.ini file