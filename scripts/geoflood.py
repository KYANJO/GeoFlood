# ----------------------------------------------
# @author:  Brian Kyanjo
# @contact: briankyanjo@u.boisestate.edu
# @date:    2022-10-16
# @version: 1.0
# @desc:    This script is to produce geflood.ini file for the geoflood model
# ------------------------------------------------

# Importing required libraries
from configparser import ConfigParser

class GeoFlooddata(object):
    # Several geoflood attributes (ignore for now)

    def __init__(self):
        self.user = {}

        self.minlevel = 0
        self.maxlevel = 0

        self.regrid_interval = 1
        self.refine_threshold = 0.5
        self.coarsen_threshold = 0.5
        self.subcycle = False
        self.output = True
        self.output_gauges = True
        self.cuda = False
        self.gravity = 9.81
        self.dry_tolerance = 1e-3
        self.earth_radius = 6371220.0
        self.coordinate_system = 1
        self.mcapa = 0
        self.verbosity = 'essential'
        self.refinement_criteria = 'value'
        self.smooth_refine = False
        self.advance_one_step = False
        self.outstyle_uses_maxlevel = True
        self.ghost_patch_pack_aux = False
        self.conservation_check = False
        self.tikz_out = False
        self.tikz_figsize = "4 4"
        self.tikz_plot_prefix = 'plot'
        self.tikz_plot_suffix = 'png'

        self.mi = 1
        self.mj = 1

        #things we didnt set in fclaw options and geoflood.py

    def write(self,rundata):
        geoflood = ConfigParser(allow_no_value=True)
        geoflood.optionxform = str # To maintain original case of comments. 

        clawdata = rundata.clawdata
        geo_data = rundata.geo_data
        refinement_data = rundata.refinement_data
        amrdata = rundata.amrdata

        user = {
        '   # User defined parameters' : None, 
        '   cuda' : self.cuda,
        '   gravity' : self.gravity,
        '   dry_tolerance' : self.dry_tolerance,
        '   earth_radius' : self.earth_radius,
        '   coordinate_system' : self.coordinate_system,
        '   mcapa' : self.mcapa
        }
        for k in self.user.keys():
            user[f'   {k:}'] = self.user[k]

        geoflood['user'] = user

        geoflood['clawpatch'] = {"\n"
        '   # Grid dimensions' : None,
        '   mx' : clawdata.num_cells[0],

        '   my' : clawdata.num_cells[1], "\n"

        '   # Number of ghost cells': None,
        '   mbc': clawdata.num_ghost, "\n"

        '   # Number of auxiliary variables' : None,
        '   maux': clawdata.num_aux,"\n"

        '   # Number of equations' : None,
        '   meqn': clawdata.num_eqn ,  "\n"   
        '   refinement-criteria' : self.refinement_criteria  
        }

        if clawdata.output_step_interval is None:
            clawdata.output_step_interval = 1


        assert self.verbosity in ['silent','essential','production','info','debug'], \
            "Error (geoflood) : Invalid value for verbosity"


        geoflood['Options'] = {"\n"

        '# Regridding options' : None,
        '   minlevel' : self.minlevel,
        '   maxlevel' : self.maxlevel,"\n"

        "   # Regrid every 'regrid_interval' time steps using threshold values":None, 
        '   regrid_interval' : self.regrid_interval,
        '   refine_threshold' : self.refine_threshold,
        '   coarsen_threshold' : self.coarsen_threshold, "\n"

        "   # Smooth refinement (around finest level)":None,
        '   smooth-refine' : self.smooth_refine,
        '   smooth-level' : self.maxlevel,"\n"


        '# Time stepping' : None, 
        '   # Final time':None,
        '   tfinal': clawdata.tfinal, "\n"     

        '   # Take a fixed time step':None,
        '   use_fixed_dt': not clawdata.dt_variable,"\n"

        "   # Initial time step for 'minlevel'":None,
        '   initial_dt': clawdata.dt_initial,"\n"

        '   # CFL constraints : Timestep will never exceed max_cfl and will ' "\n"
        '   # try not to exceed desired_cfl' : None,
        '   max_cfl': clawdata.cfl_max,
        '   desired_cfl': clawdata.cfl_desired,"\n"

        '   # 1 : Output steps  = tfinal/nout;':None,                  
        '   # 2 : not implemented;':None,                              
        '   # 3 : Take nout steps;  save files every nstep steps.':None,
        '   outstyle': clawdata.output_style,"\n" 

        '   # Used for all three out styles;  has different meaning, though':None,                                     
        '   nout': clawdata.num_output_times,"\n"
        '   # Only used if outstyle is 3':None,
        '   nstep': clawdata.output_step_interval,"\n"

        '   # Advanced time stepping' : None, 
        '   subcycle' : self.subcycle, "\n"

        '# File and console IO' : None,
        '   output' : self.output,
        '   output-gauges' : self.output_gauges,
        '   smooth-refine' : self.smooth_refine,
        '   advance-one-step' : self.advance_one_step,
        '   outstyle-uses-maxlevel' : self.outstyle_uses_maxlevel,
        '   ghost_patch_pack_aux' : self.ghost_patch_pack_aux,
        '   conservation-check' : self.conservation_check,
        '   verbosity' : self.verbosity,"\n"


        '# Domain geometry' : None,"\n"
        '   # Domain dimensions [ax,bx]x[ay,by]' : None,
        '   ax': clawdata.lower[0],
        '   bx': clawdata.upper[0],
        '   ay': clawdata.lower[1],
        '   by': clawdata.upper[1],"\n"

        '   # Block dimensions' : None,
        '   mi': self.mi,
        '   mj': self.mj,"\n"

        '   # Tikz output' : None,
        '   tikz-out': self.tikz_out,
        '   tikz-figsize': self.tikz_figsize,
        '   tikz-plot-prefix': self.tikz_plot_prefix,
        '   tikz-plot-suffix': self.tikz_plot_suffix

        }

        #mthbc
        mthbc_in  = [clawdata.bc_lower[0], clawdata.bc_upper[0], 
                     clawdata.bc_lower[1], clawdata.bc_upper[1]]
        mthbc = [0]*4
        mthbc_str = ""
        for i in range(4):
            if mthbc_in[i] in [0,'user']:        
                mthbc[i] = 0
            elif mthbc_in[i] in [1,'extrap']:    
                mthbc[i] = 1
            elif mthbc_in[i] in [2,'periodic']:  
                mthbc[i] = 2
            elif mthbc_in[i] in [3,'wall']:      
                mthbc[i] = 3
            else:
                mthbc_in[i] = mthbc[i]
            mthbc_str += " " + str(mthbc[i])

        #limiters
        lim = [0]*clawdata.num_waves
        lim_str = "" 
        for i in range(clawdata.num_waves):
            if clawdata.limiter[i] in [0,'none']:        
                lim[i] = 0
            elif clawdata.limiter[i] in [1,'minmod']:    
                lim[i] = 1
            elif clawdata.limiter[i] in [2,'superbee']:  
                lim[i] = 2
            elif clawdata.limiter[i] in [3,'vanleer']:   
                lim[i] = 3
            elif clawdata.limiter[i] in [4,'mc']:        
                lim[i] = 4
            else:
                clawdata.limiter[i] = lim[i]
            lim_str += " " + str(lim[i])

        #order 
        ord = [0]*2
        ord_str = ""
        if clawdata.order in [1,'godunov','Godunov']:
            ord[0] = 1
        if clawdata.order in [2,'Lax-Wendroff']:
            ord[0] = 2
        ord_str = str(ord[0])

        #ord_in = clawdata.transverse_waves
        if clawdata.transverse_waves in ['none' , 0]:
            ord[1] = 0
        if clawdata.transverse_waves in ['increment' , 1]:
            ord[1] = 1
        if clawdata.transverse_waves in ['all' , 2]:
            ord[1] = 2
        ord_str += " " + str(ord[1])

        # Source terms splitting:
        #   src_split == 0 or 'none'    ==> no source term (src routine never called)
        #   src_split == 1 or 'godunov' ==> Godunov (1st order) splitting used, 
        #   src_split == 2 or 'strang'  ==> Strang (2nd order) splitting used,  not recommended.

        if clawdata.source_split in [0, 'none']:
            src_split = 0
        if clawdata.source_split in [1, 'godunov', 'Godunov']:
            src_split = 1
        if clawdata.source_split in [2, 'strang', 'Strang']:
            src_split = 2

        if clawdata.output_format in [1,'ascii']:      # 'ascii' or 'binary' check line 678, "data.py", for explanation
            ascii_out = 'T'
        else:
            ascii_out = 'F'

        #print(clawdata.output_format)
        #print(ascii_out)

        geoflood['geoclaw'] = {
            '   # normal and transverse order': None,
            '   # Order of accuracy:': None,
            '   #   1 => Godunov,': None,  
            '   #   2 => Lax-Wendroff plus limiters': None,

            '   order': ord_str,"\n"

            '   # Location of capacity function in auxiliary array' : None,
            '   mcapa': clawdata.capa_index,"\n"

            '   # Source term splitting' : None,
            '   src_term': src_split,"\n"

            '   # Use an f-waves update (default : True)'
            '   use_fwaves' : clawdata.use_fwaves,"\n"

            '   # Number of waves': None,
            '   mwaves': clawdata.num_waves,"\n"


            "   # mthlim (is a vector in general, with 'mwaves' entries": None,
            '   # List of limiters to use for each wave family:': None,
            '   # Required:  len(limiter) == num_waves': None,
            '   # Some options:': None,
            "   #   0 or 'none'     ==> no limiter (Lax-Wendroff)": None,
            "   #   1 or 'minmod'   ==> minmod": None,
            "   #   2 or 'superbee' ==> superbee": None,
            "   #   3 or 'mc'       ==> MC limiter": None,
            "   #   4 or 'vanleer'  ==> van Leer": None,
            '   mthlim' : lim_str,"\n"

            '   # mthbc (=left,right,bottom,top)' : None,
            '   # Choice of BCs at xlower and xupper:':None,
            '   # 0 => user specified (must modify bcN.f to use this option)':None,
            '   # 1 => extrapolation (non-reflecting outflow)':None,
            '   # 2 => periodic (must specify this at both boundaries)':None,
            '   # 3 => solid wall for systems where q(2) is normal velocity':None,
            '   mthbc' : mthbc_str,"\n"

            '   dry_tolerance_c': geo_data.dry_tolerance,
            '   wave_tolerance_c': refinement_data.wave_tolerance, 
            '   speed_tolerance_c': refinement_data.speed_tolerance, "\n"


            '   # Output' : None,
            '   ascii-out': ascii_out
            }
            
        with open('geoflood.ini','w') as geofloodfile:
            geoflood.write(geofloodfile)


class Hydrographdata(object):

    def __init__(self):
        #  either read from file or set in setrun.py
        self.read_data = False # if true, read from file instead of setrun.py
        self.hydrograph_filename = 'filename' # name of file to read from
        self.hydrograph_type = 'discharge' # 'discharge' or 'elevation'
        self.hydrograph_variables = ['time','discharge'] # 'time','discharge','elevation'
        self.hydrograph_numrows = 0

        # hydrograph data
        self.time = []
        self.discharge = []
        self.elevation = []

        # channel data
        self.channel_width = 0.01
        self.depth = 0.1
        self.area = 0.0
        self.wetted_perimeter = 0.1
        self.friction_slope = 0.01
        self.bed_slope = 0.01
        self.froude = 1.0
        
        # initial conditions
        self.initial_velocity = 0.0
        self.initial_elevation = 0.0
        self.initial_discharge = 0.0
        self.initial_depth = 0.001

    def write(self): 
        
        # write and intial condition and hydrograph data to file
        hydrograph = open('hydrograph.data','w')
       
        # hydrograph.write('# initial condition\n')
        # hydrograph.write('# discharge (m^3/s) elevation (m) velocity (m/s)\n')
        hydrograph.write('%f %f %f %f\n' % (self.initial_discharge/max(self.channel_width,1e-8),self.initial_elevation,self.initial_velocity,self.initial_depth))

        # hydrograph.write('\n# channel width\n')
        # hydrograph.write('%f\n' % self.channel_width)

        # hydrograph.write('\n# hydrograph data\n')
        # check if hydrograph file is provided or set in setrun.py
        if self.read_data == False:
            hydrograph.write('False\n')
            if self.hydrograph_type == 'discharge':
                hydrograph.write('discharge\n')
                hydrograph.write('%d\n' % self.froude)
                hydrograph.write('%d\n' % len(self.time))
                # hydrograph.write('# time (s) discharge (m^3/s)\n')
                for i in range(len(self.time)):
                    hydrograph.write('%f %f\n' % (self.time[i],self.discharge[i]/max(self.channel_width,1e-8)))
            else:
                hydrograph.write('elevation\n')
                hydrograph.write('%d\n' % len(self.time))
                # hydrograph.write('# time (s) elevation (m)\n')                
                for i in range(len(self.time)):
                    hydrograph.write('%f %f\n' % (self.time[i],self.elevation[i]))
        else:
            hydrograph.write('True\n')
            if self.hydrograph_type == 'discharge':
                hydrograph.write('discharge\n')
                hydrograph.write('%d\n' % self.froude)
                with open(self.hydrograph_filename,'r') as hydrographfile:
                    data = [line.split() for line in hydrographfile]
                    hydrograph.write('%d\n' % len(data))
                    # hydrograph.write('# time (s) hu (m^2/s) discharge (m^3/s)\n')
                    for d in data:
                        hydrograph.write('%f %f %f\n' % (float(d[0]),float(d[1])/max(self.channel_width,1e-8),float(d[1]))) #h,hu,eta

            else:
                hydrograph.write('elevation\n')
                with open(self.hydrograph_filename,'r') as hydrographfile:
                    data = [line.split() for line in hydrographfile]
                    hydrograph.write('%d\n' % len(data))
                    # hydrograph.write('# time (s) elevation (m)\n')
                    for d in data:
                        hydrograph.write('%f %f\n' % (float(d[0]),float(d[1])))

        hydrograph.close()

# write out flowgrades data
class Flowgradesdata(object):

    def __init__(self):
        self.flowgrades = [] # list of flowgrades
        self.nflowgrades = len(self.flowgrades) # number of flowgrades

    def write(self):

        # write flowgrades data to file
        flow_grades = open('setflowgrades.data','w')
        flow_grades.write('%d\n' % len(self.flowgrades))
        for i in self.flowgrades:
            flow_grades.write(4*"%g " % tuple(i) + "\n")
        flow_grades.close()


        
