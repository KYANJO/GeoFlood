
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""

import sys
sys.path.append('../../../scripts')
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image
from geoclaw import topotools

# --------------------- Police, transformer and guage data -----------------------------------------------
malpasset_loc = "./malpasset_locs.txt"


#--------------------------
def setplot(plotdata):
#--------------------------

    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.

    """


    from visclaw import colormaps, geoplot
    from numpy import linspace

    plotdata.clearfigures()  # clear any old figures,axes,items data
    plotdata.format = "forestclaw"

    plotdata.verbose = False


    #-----------------------------------------
    # Some global kml flags
    #-----------------------------------------
    plotdata.kml_name = "Malpasset Dam"
    plotdata.kml_starttime = [1959,12,2,5,14,0]  # Date/time of event in UTC [None]
    plotdata.kml_tz_offset = 1    # Time zone offset (in hours) of event. [None]

    plotdata.kml_index_fname = "MalpassetDam"  # name for .kmz and .kml files ["_GoogleEarth"]

    # Set to a URL where KMZ file will be published.
    # plotdata.kml_publish = 'http://math.boisestate.edu/~calhoun/visclaw/GoogleEarth/kmz'

    # Add [file_name,visibility]
    # ----commented out
    # plotdata.kml_user_files.append(['malpasset_dam_validate.kml',True])

    # Cells used in setrun.py (original)
    num_cells = [54,19]
    lower = [953236.00000000, 1832407.25000000]
    upper = [959554.00000000, 1848572.75000000]

    # lower = [6.69262004,  43.40252126]
    # upper = [6.78171094,  43.55083988 ]
    #-----------------------------------------------------------
    # Figure for KML files (large view)
    # This under-resolves the finest level.
    #----------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Malpasset Dam',figno=1)
    plotfigure.show = True

    plotfigure.use_for_kml = True
    plotfigure.kml_use_for_initial_view = True

    # Latlong box used for GoogleEarth

    plotfigure.kml_xlimits = [lower[0], upper[0]]
    plotfigure.kml_ylimits = [lower[1], upper[1]]

    # Use computational coordinates for plotting
    plotfigure.kml_use_figure_limits = True
    plotfigure.kml_tile_images = False    # Tile images for faster loading.  Requires GDAL [False]


    # --------------------------------------------------
    # Resolution (should be consistent with data)
    # Refinement levels : [2,4,4,4]; max level = 5; num_cells = [55,24]
    # Aim for 1 pixel per finest level grid cell.
    # Choose a figure size that matches the coarse grid resolution.
    # Set the dpi so that figsize*dpi = finest level effective resolution.

    # If amr refinement ratios set to [0,6]; max_level = 6
    # figsize*dpi = [2,1]*16*2**6 = [2048,1024]
    mx = 16
    mi = 2
    mj = 5
    minlevel = 1
    maxlevel = 4
    p = 1
    plotfigure.kml_figsize = [2,5]  #[mx*2**p*mi,mx*2**p*mj]
    plotfigure.kml_dpi = (mi*mx*(2**maxlevel))/2  

    # --------------------------------------------------

    # Color axis : transparency below 0.1*(cmax-cmin)
    cmin = 0
    cmax = 5
    cmap = geoplot.googleearth_flooding  # transparent --> light blue --> dark blue

    # Water
    plotaxes = plotfigure.new_plotaxes('kml')
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.depth   # Plot height field h.
    # plotitem.plot_var = geoplot.surface_or_depth # mask out the lake
    plotitem.pcolor_cmap = geoplot.googleearth_flooding
    plotitem.pcolor_cmin = cmin
    plotitem.pcolor_cmax = cmax
    # plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = True      # Show patch edges

    def kml_colorbar(filename):
        geoplot.kml_build_colorbar(filename,cmap,cmin,cmax)

    plotfigure.kml_colorbar = kml_colorbar


    #-----------------------------------------------------------
    # Figure for KML files (zoomed view on region)
    #----------------------------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Transformers (zoom)',figno=2)
    plotfigure.show = True

    plotfigure.use_for_kml = True
    plotfigure.kml_use_for_initial_view = False

    # Latlong box used for GoogleEarth

    import tools
    #region_lower, region_upper, figsize = tools.region_coords(xll,xur,
    #                                                          num_cells,
    #                                                          lower,
    #                                                          upper)

    # # Get degrees per finest level :
    # dx_coarse = (upper[0]-lower[0])/num_cells[0]
    # dy_coarse = (upper[1]-lower[1])/num_cells[1]
    #
    # # Zoom region
    # mx_xlow = np.floor((xll[0] - lower[0])/dx_coarse).astype(int)
    # my_ylow = np.floor((xll[1] - lower[1])/dy_coarse).astype(int)
    # mx_zoom = np.ceil((xur[0] - xll[0])/dx_coarse).astype(int)
    # my_zoom = np.ceil((xur[1] - xll[1])/dy_coarse).astype(int)
    # xlow = lower[0] + mx_xlow*dx_coarse
    # ylow = lower[1] + my_ylow*dy_coarse

    #plotfigure.kml_xlimits = [region_lower[0],region_upper[0]]
    #plotfigure.kml_ylimits = [region_lower[1], region_upper[1]]

    # Use computational coordinates for plotting
    plotfigure.kml_use_figure_limits = True

    # --------------------------------------------------
    # plotfigure.kml_figsize = figsize*8
    # plotfigure.kml_dpi = 4*4

    # --------------------------------------------------

    plotfigure.kml_tile_images = False    # Tile images for faster loading.  Requires GDAL [False]

    # Color axis : transparency below 0.1*(cmax-cmin)
    # cmin = 0
    # cmax = 5
    # cmap = geoplot.googleearth_flooding  # transparent --> light blue --> dark blue

    # # Water
    # plotaxes = plotfigure.new_plotaxes('kml')
    # plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # plotitem.plot_var = geoplot.depth   # Plot height field h.
    # plotitem.pcolor_cmap = geoplot.googleearth_flooding
    # plotitem.pcolor_cmin = cmin
    # plotitem.pcolor_cmax = cmax

    # plot point locations 
    police, transformers, gauges, all_guages = tools.read_locations_data(malpasset_loc)
    
    # Physical model points locations (Gauges)
    gauge_points = ['P6','P7','P8','P9','P10','P11','P12','P13','P14']
    gauge_loc = []
    for i in range(len(gauges[0])):
        gauge_loc.append([gauges[1][i],gauges[2][i]])   

    # Physical model points locations (Transformers)
    transformer_points = ['A','B','c']
    transformer_loc = []
    for i in range(len(transformers[0])):
        transformer_loc.append([transformers[1][i],transformers[2][i]])

    # Physical model points locations (Police)
    police_points = ['S1','S2','S3','S4','S5','S6','S7','S8','S9','S10','S11','S12','S13','S14','S15','S16','S17']
    police_loc = []
    for i in range(len(police[0])):
        police_loc.append([police[1][i],police[2][i]])

    point_labels = transformer_points + gauge_points + police_points
    point_list = [transformer_loc, gauge_loc, police_loc]
    point_list = [item for sublist in point_list for item in sublist]
    
    # # Plot and overlay points on the plotted image
    # plt.style.use('default')
    # # im = image.imread("_plots/temp/fig1/frame0098fig1/frame0098fig1.png")
    
    # for i,pt in enumerate(point_list):
    #     lab = point_labels[i]
    #     if i>2 and i<=11:
    #         plt.plot(pt[0],pt[1],'k*',markersize=4)
    #         plt.annotate(lab, xy=pt, xytext=(pt[0]-50,pt[1]-350))
    #     elif i>11:
    #         plt.plot(pt[0],pt[1],'k*',markersize=4)
    #         plt.annotate(lab, xy=pt, xytext=(pt[0]-25,pt[1]+100))
    #     else:
    #         plt.plot(pt[0],pt[1],'r*',markersize=4)
    #         plt.annotate(lab, xy=pt, xytext=(pt[0]-5,pt[1]+100))
    
    # # plt.imshow(im)
    # plt.savefig('malpasset_points.png',dpi=300)
    
    
    #-----------------------------------------

    #-----------------------------------------
    # Figures for gauges
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Flood height', figno=300, \
                    type='each_gauge')
    plotfigure.clf_each_gauge = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'

    # Plot surface as blue curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 3  # eta?
    plotitem.plotstyle = 'b-'

    # Plot topo as green curve:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.show = True

    def gaugetopo(current_data):
        q = current_data.q
        h = q[0,:]
        eta = q[3,:]
        topo = eta - h
        return topo

    plotitem.plot_var = gaugetopo
    plotitem.plotstyle = 'g-'

    def afterframe(current_data):
        from pylab import plot, legend, xticks, floor, axis, xlabel,title
        t = current_data.t
        gaugeno = current_data.gaugeno
        # locations points
        guage_labels = police_points + transformer_points + gauge_points
        for n,pt in enumerate(guage_labels):
            if gaugeno == n:
                title(pt)

        # plot(t, 0*t, 'k')
        n = int(floor(t.max()/3600.) + 2)
        xticks([3600*i for i in range(n)], ['%i' % i for i in range(n)])
        xlabel('time (hours)')

    plotaxes.afteraxes = afterframe

    

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.parallel = False
    plotdata.print_format = 'png'           # file format
    plotdata.print_framenos = range(0,100,5)         # list of frames to print
    plotdata.print_gaugenos = 'all'         # list of gauges to print
    plotdata.print_fignos = [1,300]         # list of figures to print

    plotdata.printfigs = True              # print figures
    plotdata.overwrite = True

    plotdata.html = True                     # create html files of plots?
    plotdata.html_movie = True                    # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index

    plotdata.latex = False                    # create latex file of plots?
    #plotdata.latex_figsperline = 2           # layout of plots
    #plotdata.latex_framesperline = 1         # layout of plots
    #plotdata.latex_makepdf = False           # also run pdflatex?

    plotdata.kml = True

    return plotdata

if __name__=="__main__":
    from visclaw.plotclaw import plotclaw
    plotclaw(outdir='.',setplot=setplot,plotdir='_plots',format='forestclaw')    
